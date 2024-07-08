import math
import csv


########## FRICTIONAL RESISTANCE (R_F) CALCS ###############

def cB_calcs(displacementMass, length, beam, T):
    """
    Return Block coefficient given mass and ship dimensions
    """
    cB = displacementMass / (length * beam * T)
    if cB >= 1:
        raise Exception("Block Coefficient >= 1! Check Dimensions!")
    return displacementMass / (length * beam * T)

#block coefficient calcs
def volumetricDisplacement(length, beam, T, cB):
    """
    volumetric displacement calculation: mass / density
    """
    return length * beam * T * cB


def cP_calcs(cB, cM):
    """
    Returns prismatic coefficient, calculated from block and midship section coefficient
    """
    return cB/cM

def cF_calcs(length, v, KV):
    """
    Calculates the coefficient of friction based on length, velocity, and kinematic viscosity of water
    """
    #reynolds number calcs
    reynolds = length * v / KV #10e6 is unit converions from mm^2 to m^2
    return 0.075 / (math.log10(reynolds) - 2) ** 2

def c12_calcs(T, length): #correct
    """
    returns the c12 coefficient given the draft and length
    """
    ldRatio = T / length # draft-length ratio
    if ldRatio >= 0.05:
        return ldRatio ** 0.2228446
    if 0.02 < ldRatio and ldRatio < 0.05:
        return 48.2 * (ldRatio - 0.02) ** 2.078 + 0.479948
    if 0.02 >= ldRatio:
        return 0.479948

def c13_calcs(cStern=0):
    """
    Returns the regression-defined coefficient c13 calculated from the stern coefficient (default 0 assuming hogner stern)
    """
    return 1 + 0.003 * cStern

def lR_calcs(length, cP, lcb):
    """
    calculate lR, or run length given the length of the hull, prismatic coefficient, and longitudinal center of bouyancy
    """
    return length * ( 1 - cP + 0.06 * cP * lcb / (4 * cP - 1) )

def formFactor_calcs(c13: float, c12: float, lR: float, beam: float, cP:float, lcb: float)->float:
    """
    calculates the form factor (1 + k1)
    """
    return c13 * (0.93 + c12 * (beam / lR) ** 0.92497 * (0.95 * cP) ** -0.521448 * (1 - cP + 0.0225 * lcb) ** 0.6906)

def S_calcs(length: float, T: float, beam: float, cB: float, cM: float, cWP: float, aBT: float)->float:
    """
    Returns S, an approximation of the wetted area of the hull give length, bream, draught, and several coefficients
    """
    return length * (2 * T + beam) * math.sqrt(cM) * (0.453 + 0.4425 * cB - 0.2862 * cM - 0.003467 * beam / T + 0.3696 * cWP) + 2.38 * aBT / cB

def rF_calcs(rho, cF, S, v)->float:
    """
    computes the total resistance due to friction on the ship hull
    """
    return 0.5 * rho * cF * S * (v ** 2)



########## APPENDANGE RESISTANCE (R_APP) CALCS ###############


def rAPP_calcs(rho, v, sAPP, cF, flowAppendage = 1.5)->float:
    """
    calculates resistance from appendages. Default of appendage flow is 1.5
    """
    return 0.5 * rho * (v ** 2) * sAPP * flowAppendage * cF
#resistance from appendages. 1 is the water density


########## Wave RESISTANCE (R_W) CALCS ###############

def c7_calcs(beam, length)->float:
    """
    Calculates the regression-computed coefficient c7, derived from beam and length
    """
    if (beam / length) < 0.11:
        return 0.229577 * (beam / length) ** (1/3)
    elif 0.11 <= beam / length and beam / length <= 0.25:
        return beam / length
    elif 0.25 < beam / length:
        return 0.5 - 0.0625 * length / beam

def c3_calcs(aBT, beam, TF, T, hB)->float:
    """
    Calculates the regression-computed coefficient c3, which models the influence of a bulbous bow on wave resistance
    """
    return 0.56 * aBT ** 1.5 / (beam * T * (0.31 * math.sqrt(aBT) + TF - hB))

def c2_calcs(c3, aBT)->float:
#c2 is wave resistance coefficient due to bulbous bow. defaults to 1 for non-bulbous bows
    if aBT == 0: return 1
    else: return math.exp(-1.89 * math.sqrt(c3))

def c5_calcs(aT, beam, T, cM)->float:
    """
    Calculates c5, the wave resistance coefficient due to transom sterns. defaults to 1 for non-transom sterns
    """
    if aT == 0: return 1
    else: return 1 - 0.8 * aT / (beam * T * cM)

def nabla_calcs(cB, length, beam, T)->float:
    """
    Calculates the moulded displacement volume
    """
    return cB  * length * beam * T

def iE_calcs(length, beam, cWP, cP, lcb, lR, nabla)->float:
    """
    calculates the value iE, the half angle of entrance (angle of the waterline at the bow with reference to the center plane)
    This value can be manually inputted based on the actual design, but estimated via regression analysis for user simplicity
    """
    return 1 + 89 * math.exp(-1 * (length / beam) ** 0.80856 * (1 - cWP) ** 0.30484 * (1 - cP - 0.0225 * lcb) ** 0.6367 
    * (lR / beam) ** 0.34574 * (100 * nabla / length ** 3) ** 0.16302)

def c1_calcs(c7, T, beam, iE)->float:
    """
    Calculates the regression-derived coefficient c1, based on c7, the draught and beam, and the half-angle of entrance
    """
    return 2223105 * c7 ** 3.78613 * (T / beam) ** 1.07961 * (90 - iE) ** -1.37565

def froude_length_calcs(v, length, G)->float:
    """
    calculates the froude number given the velocity and waterline length
    """
    return v / math.sqrt(length * G)

def lamdba_calcs(length, beam, cP)->float:
    """
    calculates lambda, a regression-derived coefficient relating to wave resistance
    """
    if (length / beam) < 12: return 1.446 * cP - 0.03 * length / beam
    else: return 1.446 * cP - 0.36

def c16_calcs(cP)->float:
    """
    Calculates c16, a regression-derived coefficient based on the prismatic coefficient
    """
    if cP <= 0.80: return 8.07981 * cP - 13.8763 * cP ** 2 + 6.984388 * cP ** 3
    else: return 1.73014 - 0.7067

def m1_calcs(length, T, nabla, beam, c16)->float:
    """
    calculates m1, a regression-derived parameter relating to wave resistances
    """
    return 0.0140407 * length / T - 1.75254 * nabla ** (1/3) / length - 4.79323 * beam / length - c16

def c15_calcs(length, nabla)->float:
    """
    calculates c15, a regression-derived coefficient based on the ratio of the length cubed and nabla
    """
    length_nabla = (length ** 3) / nabla
    if length_nabla < 512: return -1.69385
    elif length_nabla > 1727: return 0
    else: return -1.69385 + (length / nabla ** (1.3) - 8.0)/2.36

def m2_calcs(c15, cP, Fn)->float:
    """
    calculates m2, a regression-derived parameter relating to wave resistances
    """
    return c15 * cP ** 2 * math.exp(-0.1 * Fn ** -2)

def rW_calcs(c1, c2, c5, nabla, rho, G, m1, Fn, m2, lambda_w)->float:
    """
    calculates wave resistance based on previously calculated parameters and coefficients 
    """
    d = -0.9
    return c1 * c2 * c5 * nabla * rho * G * math.exp(m1 * Fn ** d + m2 * math.cos(lambda_w * Fn ** -2))

########## BULBOUS BOW RESISTANCE (R_B) CALCS ###############

def rB_calcs(aBT, hB, G, TF, v, rho)->float:
    """
    calculates the resistance of a bulbous bow given the ship speed and dimensions of the bow
    """
    #measure of emergence of the bow
    p_B = 0.56 * math.sqrt(aBT) / (TF - 1.5 * hB)

    #Froude Number based on immersion of the bow
    Fni = v / math.sqrt(G * (TF - hB - 0.25 * math.sqrt(aBT)) + 0.15 * v ** 2)

    #0 if no bulbous bow 
    if aBT == 0: return 0
    else: return 0.11 * math.exp(-3 * p_B ** -2) * Fni ** 3 * aBT ** 1.5 * rho * G / (1 + Fni ** 2)



########## TRANSOM STERN RESISTANCE (R_TR) CALCS ###############

def rTR_calcs(aT, v, G, beam, cWP, rho)->float:
    """
    calculates the resistance caused by a tramson stern based on the dimensions of the ship and the stern
    """
    #only run calc if transom stern exists
    if aT > 0:
        #Froude Number based on transom immersion
        FnT = v / math.sqrt(2 * G * aT / (beam + beam * cWP))
        if FnT < 5: c6 = 0.2 * (1 - 0.2 * FnT)
        else: c6 = 0
        return 0.5 * rho * v ** 2 * aT * c6
    else: return 0



########## MODEL-SHIP CORRELATION RESISTANCE (R_A) CALCS ###############

def cA_calcs(TF, length, cB, c2)->float:
    """
    calculates the correlation allowance coefficient
    """
    #c4 factors in the tramsom stern
    if TF / length <= 0.04: c4 = TF / length
    else: c4 = 0.04
    #
    return 0.006 * (length + 100) ** -0.16 - 0.00205 + 0.003 * math.sqrt(length / 7.5) * cB ** 4 * c2 * (0.04 - c4)


def rA_calcs(v, cA, rho, S)->float:
    """ 
    Calculates the model-ship correlation resistance
    """
    return 0.5 * rho * v **2 * S * cA


########## TOTAL RESISTANCE CALCS ###############


def rTotal_calcs(rF, formFactor, rAPP, rW, rB, rTR, rA)->float:
    """
    calculates total resistance force on the ship (as newtons) using previously calculated factors
    rTotal = cF * (formFactor) + rAPP + rW + rB + rTR + rA
    """
    #total resistance (newtons)
    return rF * (formFactor) + rAPP + rW + rB + rTR + rA

def externalPower_calcs(rTotal, v)->float:
    """
    calculates the external/effective power applied to the ship as a whole, NOT the shaft power, which is the power applied to the shafts
    based on the formula that power equals force * velocity
    """
    return rTotal * v


def cV_calcs(formFactor, cF, cA)->float:
    """
    Calculates the viscious coefficient, which factors in the friction of the water on the ship as well as
    the influence of hull form on viscous pressure drag.
    """
    return formFactor * cF + cA




########## PROPELLER-VARIABLE CALCS (w, t, eta_R) ###############
#any calculations that differ depending on the number of propellors for a ship. 
#Holtrop & Mennen only account for 1 or 2 propellers, 
#so anything after that is assumed to have similar efficiency to a twin-screw arrangement
#also assume aft draught (T_A)is the same as the amidships draught (T)

def c8_calcs(beam, T, S, length, dProp)->float:
    """
    calculate the propeller-variable coefficient c8
    """
    if beam / T < 5:
        return beam * S / (length * dProp * T)
    else:
        return S * (7 * beam / T - 25) / (length * dProp * (beam / T - 3))



def c9_calcs(c8)->float:
    """
    calculate the propeller-variable coefficient c9
    """
    if c8 < 28: return c8
    else:
        return 32 - 16 / (c8 - 24)



def c10_calcs(beam, length)->float:
    """
    calculate the regression-derived propeller-variable coefficient c10
    """
    if length / beam > 5.2:
        return beam / length
    else:
        return 0.25 - 0.003328402 * (beam / length - 0.134615385)



def c11_calcs(T, dProp)->float:
    """
    calculate the regression-derived propeller-variable coefficient c11 based on the draft (at the rear of the ship) and the propeller diameter
    """
    if T/dProp < 2:
        return T / dProp
    else:
        return 0.0833333 * (T / dProp) ** 3 + 1.33333



def cP1_calcs(cP, lcb)->float:
    """
    Calculate the regression-derived coefficient cP1, a prismatic coefficient variable accounting for longitudinal centre of bouyancy
    """
    return 1.45 * cP - 0.315 - 0.0225 * lcb



def c20_calcs(cStern)->float:
    """
    Calculate the regression-derived coefficient c20 based on the stern coefficient
    """
    return 1.0 + 0.015 * cStern



def c19_calcs(cB, cP, cM)->float:
    """
    Calculate the regression-derived coefficient c19 based on the hull form coefficients of the ship
    """
    if cP < 0.7:
        return 0.12997 / (0.95 - cB) - 0.11056 / (0.95 - cP)
    else:
        return 0.18567 / (1.3571 - cM) - 0.71276 + 0.38648 * cP



def wt_calcs(c9, c20, cV, length, T, c11, cP1, beam, c19, dProp, cStern, cB, numPropellers, c10)->float:
    """
    Calculates the wake fraction (w) and thrust deduction coefficients (t).
    """
    if numPropellers == 1:
        #assume conventional stern (not an open stern)
        #not using revised method

        #use revision for calculating single screw conventional stern, Holtrop 1984
        w = c9 * c20 * cV * (length / T) * (0.050776 + 0.93405 * c11 * (cV / (1 - cP1))) + (0.27915 * c20 * math.sqrt(beam / (length * (1 - cP1)))) + c19 * c20
        #still use Holtrop 1978 calculation
        td = 0.001979 * length / (beam - beam * cP1) + 1.0585 * c10 
        - 0.00524 - 0.0148 * dProp ** 2 / (beam * T) + 0.0015 * cStern
    else: #2 propellers (extrapolated to 2+ propellers)
        w = 0.3095 * cB + 10 * cV * cB - 0.23 * dProp / math.sqrt(beam * T)
        td = 0.325 * cB - 0.1885 * dProp / math.sqrt(beam * T)
    return w, td




########## BLADE AREA RATIO (A_E/A_O) CALCS ###############

def K_calcs(numPropellers)->float:
    """
    K is a constant based on number of propellers. Only determined for single or double screw ships; unreliable for more than 2 screws

    """
    if numPropellers == 1: return 0.2
    elif numPropellers > 1: return 0.1
    else: raise Exception("Ship cannot have less than 1 propeller!")



def hShaft_calcs(T, propKeelClearance, dProp)->float:
    """
    returns the depth of the shaft centerline, calculated as draught minus prop radius + prop keel clearance (at the stern)
    """
    return T - (propKeelClearance + dProp / 2)





def pitch_calcs(v, w, n)->float:
    """
    calculates propeller pitch, or how far the propeller goes in one revolution
    formula: pitch velocity * wake field coefficient (to account for propeller slip) / rotational frequency
    """
    return v * (1 - w) / n

def propThrust_calcs(rTotal, td)->float:
    """
    returns propeller thrust the thrust generated by the propellers
    effective thrust = total resistance
    effective thrust = prop thrust * thrust deduction coefficienct (1 - t)
    prop thrust = total resistance / thrust deduction coefficient
    """
    return rTotal / (1 - td)



def bladeAreaRatio_calcs(K, numBlades, propThrust, dProp, rho, G, hShaft)->float:
    """
    Calculates the ratio of the area of the blades 
    original equation: bladeAreaRatio = K + (1.3 + 0.3 * Z) * thrust / (dProp ** 2 * (p_o + rho * G * h - p_v))
    modified to replace p_o - p_v with 99047 N/m^2, which holds true for seawater at 15 degrees celsius

    """
    return  K + (1.3 + 0.3 * numBlades) * propThrust / (dProp ** 2 * (99047 + rho * G * hShaft))


def etaR_calcs(numPropellers, bladeAreaRatio, cP, lcb, pitch, dProp)->float:
    #calculate eta_R here since it is dependent on A_E/A_O for the case of a single screw ship
    if numPropellers == 1:
        return 0.9922 - 0.05908 * bladeAreaRatio + 0.07424 * (cP - 0.0225 * lcb)
    else:
        return 0.9737 + 0.111 * (cP - 0.0225 * lcb) + 0.06325 * pitch / dProp

########## OPEN WATER PROPELLER EFFICIENCY CALCS (eta_o) ###############

def chordLength_calcs(bladeAreaRatio, dProp, numBlades)->float:
    """
    calculate c_0.75, or the chord length of the propeller
    """
    return 2.073 * (bladeAreaRatio) * dProp / numBlades

def tc_calcs(numBlades, c_075)->float:
    """
    returns the thickness-chord length ratio of the propeller
    """
    return (0.0185 - 0.00125 * numBlades) / c_075

def deltaCD_calcs(tc_075, c_075, k_p)->float:
    """
    Calculates the different in drag coefficients of the hull profile section
    """
    return (2 + 4 * tc_075) * (0.003605 - (1.89 + 1.62 * math.log(c_075 / k_p)) ** -2.5)

def ktb_calcs(propThrust, rho, dProp, n)->float:
    """
    Calculates the thrust coefficient of a model (K_T_B or K_TM) based on the ITTC 1978 specification formula.


    """

    #ITTC 1978 formulas
    return propThrust / (rho * dProp ** 4 * n ** 2)




def kt_calcs(K_T_B, delta_CD, pitch, c_075, numBlades, dProp)->float:
    """
    calculates the ship thrust coefficient using the Holtrop-Mennen formulation
    """
    return K_T_B + delta_CD * 0.3 * (pitch * c_075 * numBlades) / dProp ** 2


def j_calcs(v, w, td, n, dProp)->float:
    """
    calculate the advance ratio of the propeller J
    """
    return v * (1 - w * td) / (n * dProp)


#can't figure out torque
#K_Q_B = propTorque / (rho * dProp ** 5 * n ** 2)
#K_Q = K_Q_B - delta_CD * 0.25 * (c_075 * numBlades) / dProp
#give up, can't figure way to derive torque
#eta_o = J * K_T /(2 * math.pi * K_Q)

def cTH_calcs(K_T, J)->float:
    """
    calculates the thrust coefficient, using the formulas provided by ITTC 1978 proceedings
    """
    return (K_T / J ** 2) * 8 / math.pi


def etao_calcs(cTH, trueEfficiencyCoefficient = 0.7)->float:
    """
    calculates eta_o, the propeller efficiency. Use a textbook ideal propeller efficiency formula, which does overestimate propeller efficincy
    0.7 coefficient derived from fitting to the example ship provided by Holtrop & Mennen
    """
    #ideal propeller efficiency / ideal eta_o calculation, is a drastic overestimation of eta_o:
    eta_o = 2 / (1 + math.sqrt(1 + cTH))

    #multiply ideal propeller efficiency by 0.7 to reflect real life conditions
    #still results in optimistic measurements
    #eta_o *= trueEfficiencyCoefficient
    return eta_o
########## TOTAL SHAFT POWER CALCS ###############

def shaftPowerCalcs(P_E, eta_R, eta_o, eta_S, td, w)->float:
    """
    Calculates the required shaft power given the external power and all the efficiency coefficients
    """
    return P_E / (eta_R * eta_o * eta_S * (1 - td)/(1 - w))


def HoltropMennenPowerCalculation(length, beam, T, displacementMass, v,
                                  lcb = 0, cM = 0.95, sAPP = 0, cWP = 0.7, aBT = 0,
                                  hB = 4, aT = 0, numPropellers = 2, dProp = 3.5,
                                  numBlades = 3, n = 3, propKeelClearance = 0.2, trueEfficiencyCoefficient = 0.7):
    """
    Full calculation of resistance and shaft power. length, beam, draft, displacement, and velocity are required, all other params are optional.
    All units are in metric
    PARAMETERS:
    length: waterline length of ship in meters
    beam: beam of the ship (at its widest) in meters
    T: the average draft/depth of the ship at the given displacement in meters. program assumes uniform draft.
    displacementMass: the displacement of the ship at a given load in metric tons. Trial measurements are ideal, but other displacements are fine.
    v: the speed of the ship, in meters per second
    lcb: longitudinal center of bouynacy (% of ship's length in front of amidships [0.5 L]). Default is 0, or perfectly amidships
    cM: midship section coefficient (midship section area / beam * draft). Merchant ships are 0.9, Bismarck was 0.97
        Tegethoff (1876) was 0.82, Minas Geraes was 0.967, high speed destroyer should have 0.8. Default is 0.95
    sAPP: wetted area appendage. Appendages are any underwater structures protruding from the hull, like the rudder and propellors. Default is 0 if unsure.
    cWP: water plane coefficient (waterplane area / beam * length). cWP of 1 is a rectangle. Larger ships have a slightly lower cWP than smaller ones.
        Farraguts were 0.744, post-WWII destroyers were 0.68, Bismarck was 0.66, Minas Geraes was 7.12. Default of 0.7.
    aBT: cross-sectional area of the bulbous bow in square meters. 0 for non bulbous bow (default)
    hB: height of the center of the bulbous bow above keel line in meters. default of 4, not used if aBT = 0
    aT: immersed area of the transverse area of the transom at zero speed in square meters. Default is 0 for no tramson stern
    numPropellers: number of shafts driving propellers (thrusters exclusively for steering don't count).
    dProp: average diameter of the propellers (meters). Default of 3.5 meters
    numBlades: number of blades on each propeller. Default of 3
    n: the average shaft speed in rotations per second or hertz. Default of 3 rotations per second. 
        HMS Hood was 210 (3.5 s^-1), Bismarck was 265 rpm (4.416 s^-1), USS Massachussets was 185 rpm (3.083 rps),
        RMS Carmania was 185 rpm, Minas Geraes was 147rpm (2.45 rps)
    propKeelClearance: how many meters between the tip of a propellor at its lowest and the keel line. Default to 0.2
    """
    ####### INPUT VALUES ############


    #TF is the forward draft. Ships are assumed to be large and slow enough that the difference between the amidships and forward draft is negligible
    TF = T
    ########## CONSTANTS ###############

    G = 9.81 #acceleration due to gravity (m/s^2)
    rho  = 1025  #density of fluid (salt water in this case) kg / m^3
    KV = 11.8987e-7 #kinematic viscoscity (m^2/s). Default value of 11.8987e-7 for water at 16 celsius

    ########## APPROXIMATED PARAMETERS ###############
    k_p = 0.00003 #propeller roughness (for 1980s props, may be higher with fouled or older propellers, but assume 0.00003 for simplicity)
    #eta_s is the sea efficincy coefficient, which accounts for sea currents, water quality, and fouling
    eta_S = 0.99 #coefficient ideal conditions (completely calm water, 15 degrees saltwater, clean hull + propeller)
    cCrossSection = 0.9 #cross-section coefficient. 1 is a rectangle, 0.5 is a triangle. 0.9 default assuming u-shaped hull
    flowAppendage = 1.5 #weighted average approximation of the appendage flow coefficient. See (Holthrop and Mennen 167) for a more accurate calculation

    #RESISTANCE CALCULATIONS
    #wetted area of the hull
    cB = cB_calcs(displacementMass, length, beam, T)
    cP = cP_calcs(cB, cM) #prismatic coefficient calculation
    cF = cF_calcs(length, v, KV) #cF is coefficient of friction
    cStern = 0 # 10 for U-shaped section with hogner stern, 0 for normal section, -10 for v-shaped sections
    c12 = c12_calcs(T, length)
    c13 = c13_calcs(cStern)
    lR = lR_calcs(length, cP, lcb)
    formFactor = formFactor_calcs(c13, c12, lR, beam, cP, lcb)
    S = S_calcs(length, T, beam, cB, cM, cWP, aBT)
    rF = rF_calcs(rho, cF, S, v)
    rAPP = rAPP_calcs(rho, v, sAPP, cF, flowAppendage)
    c7 = c7_calcs(beam, length)
    #influence of bulbous bow on wave resistance
    c3 = c3_calcs(aBT, beam, TF, T, hB)
    c2 = c2_calcs(c3, aBT)
    c5 = c5_calcs(aT, beam, T, cM)
    #moulded displacement volume (nabla)
    nabla = nabla_calcs(cB, length, beam, T)
    iE = iE_calcs(length, beam, cWP, cP, lcb, lR, nabla)
    c1 = c1_calcs(c7, T, beam, iE)
    #Froude Number, based on waterline length
    Fn = froude_length_calcs(v, length, G)
    print(f"Froude Number: {Fn}")
    lambda_w = lamdba_calcs(length, beam, cP)
    c16 = c16_calcs(cP)
    m1 = m1_calcs(length, T, nabla, beam, c16)
    c15 = c15_calcs(length, nabla)
    m2 = m2_calcs(c15, cP, Fn)
    rW = rW_calcs(c1, c2, c5, nabla, rho, G, m1, Fn, m2, lambda_w)
    rB = rB_calcs(aBT, hB, G, TF, v, rho)
    rTR = rTR_calcs(aT, v, G, beam, cWP, rho)
    cA = cA_calcs(TF, length, cB, c2)
    rA = rA_calcs(v, cA, rho, S)
    rTotal = rTotal_calcs(rF, formFactor, rAPP, rW, rB, rTR, rA)
    #total resistance deviates 4% from Holthrop + Mennen's example (underestimation of resistance)
    P_E = externalPower_calcs(rTotal, v)
    #power deviates ~2% from Holthrop + Mennen's example (underestimation)
    cV = cV_calcs(formFactor, cF, cA)
    #SHAFT POWER CALCULATIONS
    c8 = c8_calcs(beam, T, S, length, dProp)
    c9 = c9_calcs(c8)
    c10 = c10_calcs(beam, length)
    c11 = c11_calcs(T, dProp)
    cP1 = cP1_calcs(cP, lcb)
    c20 = c20_calcs(cStern)
    c19 = c19_calcs(cB, cP, cM)
    w, td = wt_calcs(c9, c20, cV, length, T, c11, cP1, beam, c19, dProp, cStern, cB, numPropellers, c10)
    K = K_calcs(numPropellers)
    hShaft = hShaft_calcs(T, propKeelClearance, dProp)
    pitch = pitch_calcs(v, w, n)
    propThrust = propThrust_calcs(rTotal, td)
    bladeAreaRatio = bladeAreaRatio_calcs(K, numBlades, propThrust, dProp, rho, G, hShaft)
    eta_R = etaR_calcs(numPropellers, bladeAreaRatio, cP, lcb, pitch, dProp)
    #chord length
    c_075 = chordLength_calcs(bladeAreaRatio, dProp, numBlades)
    tc_075 = tc_calcs(numBlades, c_075)

    #difference in drag coefficients of the profile section
    delta_CD = deltaCD_calcs(tc_075, c_075, k_p)
    K_T_B = ktb_calcs(propThrust, rho, dProp, n)
    K_T = kt_calcs(K_T_B, delta_CD, pitch, c_075, numBlades, dProp)
    J = j_calcs(v, w, td, n, dProp)
    cTH = cTH_calcs(K_T, J)
    eta_o = etao_calcs(cTH, trueEfficiencyCoefficient = 0.7)
    shaftPower = shaftPowerCalcs(P_E, eta_R, eta_o, eta_S, td, w)
    shaftPower *= 1 / trueEfficiencyCoefficient
    return shaftPower

def main():
    with open("HoltropMennenTest.csv", 'r') as csvfile:
        # creating a csv reader object
        csvreader = csv.reader(csvfile)
        fields = next(csvreader)
        for row in csvreader:
            for idx in range(1, 14): 
                if row[idx] != '': row[idx] = float(row[idx].replace(',', ''))
            name, length, beam, draft, displacement, speed, numShafts, numBlades = row[0], row[2], row[3], row[4], row[5], row[9], row[10], row[11]
            speed /= 1.944
            if row[7] != '': cM = row[7]
            else: cM = 0.95
            if row[8] != '': cWP = row[8]
            else: cWP = 0.7
            if row[12] != '': dProp = row[12]
            else: dProp = 3.5
            if row[13] != '': propSpeed = row[13]
            else: propSpeed = 3
            shaftPower = HoltropMennenPowerCalculation(length, beam, draft, displacement, speed, cM = cM, cWP = cWP,
                                numPropellers = numShafts, dProp = dProp,
                                numBlades = numBlades, n = propSpeed)
            print(f"{name}: {round(shaftPower/1000)} kW")
            


if __name__ == "__main__":
    main()

#sources: overall calculation (Holthrop and Mennen 1978): https://repository.tudelft.nl/islandora/object/uuid:ee370fed-4b4f-4a70-af77-e14c3e692fd4/datastream/OBJ/download
#coefficient of friction + reynolds number calculation: https://repository.tudelft.nl/islandora/object/uuid%3A16d77473-7043-4099-a8c6-bf58f555e2e7
#prismatic coefficient (cP): https://www.nautilusshipping.com/form-coefficient-of-ship
#destroyer ship coefficient numbers: https://www.jstor.org/stable/pdf/44893811.pdf
#calculating open water propeller efficiency (eta_o): ITTC 1978 https://ittc.info/media/1834/75-02-03-014.pdf
#open water propeller efficiency: https://apps.dtic.mil/sti/tr/pdf/ADA132414.pdf
#open water propeller efficiency: https://www.man-es.com/docs/default-source/document-sync/basic-principles-of-ship-propulsion-eng.pdf
#Wagignen B-Series Polynomial