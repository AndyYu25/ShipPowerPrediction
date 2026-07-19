import math
import csv


########## FRICTIONAL RESISTANCE (R_F) CALCS ###############

def cB_calcs(displacementMass: float, length: float, beam: float, T: float) -> float:
    """Return Block coefficient given mass and ship dimensions.

    :param displacementMass: Displacement mass of the ship.
    :type displacementMass: float
    :param length: Waterline length of the ship.
    :type length: float
    :param beam: Beam of the ship (at its widest point).
    :type beam: float
    :param T: Draft/depth of the ship.
    :type T: float
    :returns: Block coefficient of the ship.
    :rtype: float
    :raises Exception: If the computed block coefficient is >= 1, indicating invalid dimensions.
    """
    cB = displacementMass / (length * beam * T)
    if cB >= 1:
        raise Exception("Block Coefficient >= 1! Check Dimensions!")
    return displacementMass / (length * beam * T)

#block coefficient calcs
def volumetricDisplacement(length: float, beam: float, T: float, cB: float) -> float:
    """Volumetric displacement calculation: mass / density.

    :param length: Waterline length of the ship.
    :type length: float
    :param beam: Beam of the ship.
    :type beam: float
    :param T: Draft of the ship.
    :type T: float
    :param cB: Block coefficient of the ship.
    :type cB: float
    :returns: Volumetric displacement of the hull.
    :rtype: float
    """
    return length * beam * T * cB


def cP_calcs(cB: float, cM: float) -> float:
    """Return prismatic coefficient, calculated from block and midship section coefficient.

    :param cB: Block coefficient.
    :type cB: float
    :param cM: Midship section coefficient.
    :type cM: float
    :returns: Prismatic coefficient.
    :rtype: float
    """
    return cB/cM

def cF_calcs(length: float, v: float, KV: float = 11.8987e-7) -> float:
    """Calculate the coefficient of friction based on length, velocity, and kinematic viscosity of water.

    :param length: Waterline length of the ship.
    :type length: float
    :param v: Ship speed.
    :type v: float
    :param KV: Kinematic viscosity of water, defaults to 11.8987e-7.
    :type KV: float, optional
    :returns: Coefficient of friction.
    :rtype: float
    """
    #reynolds number calcs
    reynolds = length * v / KV #10e6 is unit converions from mm^2 to m^2
    return 0.075 / (math.log10(reynolds) - 2) ** 2

def c12_calcs(T: float, length: float) -> float:
    """Return the c12 coefficient given the draft and length.

    :param T: Draft of the ship.
    :type T: float
    :param length: Waterline length of the ship.
    :type length: float
    :returns: The c12 coefficient.
    :rtype: float
    """
    ldRatio = T / length # draft-length ratio
    if ldRatio >= 0.05:
        return ldRatio ** 0.2228446
    if 0.02 < ldRatio and ldRatio < 0.05:
        return 48.2 * (ldRatio - 0.02) ** 2.078 + 0.479948
    if 0.02 >= ldRatio:
        return 0.479948

def c13_calcs(cStern: float = 0) -> float:
    """Return the regression-defined coefficient c13 calculated from the stern coefficient.

    :param cStern: Stern shape coefficient, defaults to 0 (assumes hogner stern).
    :type cStern: float, optional
    :returns: The c13 coefficient.
    :rtype: float
    """
    return 1 + 0.003 * cStern

def lR_calcs(length: float, cP: float, lcb: float) -> float:
    """Calculate lR, or run length, given the length of the hull, prismatic coefficient, and longitudinal center of bouyancy.

    :param length: Waterline length of the ship.
    :type length: float
    :param cP: Prismatic coefficient.
    :type cP: float
    :param lcb: Longitudinal center of bouyancy.
    :type lcb: float
    :returns: Run length (lR).
    :rtype: float
    """
    return length * ( 1 - cP + 0.06 * cP * lcb / (4 * cP - 1) )

def formFactor_calcs(c13: float, c12: float, lR: float, beam: float, cP:float, lcb: float)->float:
    """Calculate the form factor (1 + k1).

    :param c13: The c13 coefficient.
    :type c13: float
    :param c12: The c12 coefficient.
    :type c12: float
    :param lR: Run length.
    :type lR: float
    :param beam: Beam of the ship.
    :type beam: float
    :param cP: Prismatic coefficient.
    :type cP: float
    :param lcb: Longitudinal center of bouyancy.
    :type lcb: float
    :returns: The form factor (1 + k1).
    :rtype: float
    """
    return c13 * (0.93 + c12 * (beam / lR) ** 0.92497 * (0.95 * cP) ** -0.521448 * (1 - cP + 0.0225 * lcb) ** 0.6906)

def S_calcs(length: float, T: float, beam: float, cB: float, cM: float, cWP: float, aBT: float)->float:
    """Return S, an approximation of the wetted area of the hull given length, beam, draught, and several coefficients.

    :param length: Waterline length of the ship.
    :type length: float
    :param T: Draft of the ship.
    :type T: float
    :param beam: Beam of the ship.
    :type beam: float
    :param cB: Block coefficient.
    :type cB: float
    :param cM: Midship section coefficient.
    :type cM: float
    :param cWP: Waterplane coefficient.
    :type cWP: float
    :param aBT: Cross-sectional area of the bulbous bow.
    :type aBT: float
    :returns: Approximated wetted area of the hull.
    :rtype: float
    """
    return length * (2 * T + beam) * math.sqrt(cM) * (0.453 + 0.4425 * cB - 0.2862 * cM - 0.003467 * beam / T + 0.3696 * cWP) + 2.38 * aBT / cB

def rF_calcs(rho: float, cF: float, S: float, v: float) -> float:
    """Compute the total resistance due to friction on the ship hull.

    :param rho: Density of the fluid.
    :type rho: float
    :param cF: Coefficient of friction.
    :type cF: float
    :param S: Wetted area of the hull.
    :type S: float
    :param v: Ship speed.
    :type v: float
    :returns: Frictional resistance.
    :rtype: float
    """
    return 0.5 * rho * cF * S * (v ** 2)



########## APPENDANGE RESISTANCE (R_APP) CALCS ###############


def rAPP_calcs(rho: float, v: float, sAPP: float, cF: float, flowAppendage: float = 1.5) -> float:
    """Calculate resistance from appendages.

    :param rho: Density of the fluid.
    :type rho: float
    :param v: Ship speed.
    :type v: float
    :param sAPP: Wetted area of the appendages.
    :type sAPP: float
    :param cF: Coefficient of friction.
    :type cF: float
    :param flowAppendage: Appendage flow coefficient, defaults to 1.5.
    :type flowAppendage: float, optional
    :returns: Appendage resistance.
    :rtype: float
    """
    return 0.5 * rho * (v ** 2) * sAPP * flowAppendage * cF
#resistance from appendages. 1 is the water density


########## Wave RESISTANCE (R_W) CALCS ###############

def c7_calcs(beam: float, length: float) -> float:
    """Calculate the regression-computed coefficient c7, derived from beam and length.

    :param beam: Beam of the ship.
    :type beam: float
    :param length: Waterline length of the ship.
    :type length: float
    :returns: The c7 coefficient.
    :rtype: float
    """
    if (beam / length) < 0.11:
        return 0.229577 * (beam / length) ** (1/3)
    elif 0.11 <= beam / length and beam / length <= 0.25:
        return beam / length
    elif 0.25 < beam / length:
        return 0.5 - 0.0625 * length / beam

def c3_calcs(aBT: float, beam: float, TF: float, T: float, hB: float) -> float:
    """Calculate the regression-computed coefficient c3, which models the influence of a bulbous bow on wave resistance.

    :param aBT: Cross-sectional area of the bulbous bow.
    :type aBT: float
    :param beam: Beam of the ship.
    :type beam: float
    :param TF: Forward draft of the ship.
    :type TF: float
    :param T: Draft of the ship.
    :type T: float
    :param hB: Height of the center of the bulbous bow above the keel line.
    :type hB: float
    :returns: The c3 coefficient.
    :rtype: float
    """
    return 0.56 * aBT ** 1.5 / (beam * T * (0.31 * math.sqrt(aBT) + TF - hB))

def c2_calcs(c3: float, aBT: float) -> float:
#c2 is wave resistance coefficient due to bulbous bow. defaults to 1 for non-bulbous bows
    if aBT == 0: return 1
    else: return math.exp(-1.89 * math.sqrt(c3))

def c5_calcs(aT: float, beam: float, T: float, cM: float) -> float:
    """Calculate c5, the wave resistance coefficient due to transom sterns.

    :param aT: Immersed area of the transverse transom.
    :type aT: float
    :param beam: Beam of the ship.
    :type beam: float
    :param T: Draft of the ship.
    :type T: float
    :param cM: Midship section coefficient.
    :type cM: float
    :returns: The c5 coefficient, defaults to 1 for non-transom sterns.
    :rtype: float
    """
    if aT == 0: return 1
    else: return 1 - 0.8 * aT / (beam * T * cM)

def nabla_calcs(cB: float, length: float, beam: float, T: float) -> float:
    """Calculate the moulded displacement volume.

    :param cB: Block coefficient.
    :type cB: float
    :param length: Waterline length of the ship.
    :type length: float
    :param beam: Beam of the ship.
    :type beam: float
    :param T: Draft of the ship.
    :type T: float
    :returns: Moulded displacement volume (nabla).
    :rtype: float
    """
    return cB  * length * beam * T

def iE_calcs(length: float, beam: float, cWP: float, cP: float, lcb: float, lR: float, nabla: float) -> float:
    """Calculate the value iE, the half angle of entrance.

    This value can be manually inputted based on the actual design, but is estimated via
    regression analysis here for user simplicity.

    :param length: Waterline length of the ship.
    :type length: float
    :param beam: Beam of the ship.
    :type beam: float
    :param cWP: Waterplane coefficient.
    :type cWP: float
    :param cP: Prismatic coefficient.
    :type cP: float
    :param lcb: Longitudinal center of bouyancy.
    :type lcb: float
    :param lR: Run length.
    :type lR: float
    :param nabla: Moulded displacement volume.
    :type nabla: float
    :returns: Half angle of entrance (iE), the angle of the waterline at the bow with reference to the center plane.
    :rtype: float
    """
    return 1 + 89 * math.exp(-1 * (length / beam) ** 0.80856 * (1 - cWP) ** 0.30484 * (1 - cP - 0.0225 * lcb) ** 0.6367 
    * (lR / beam) ** 0.34574 * (100 * nabla / length ** 3) ** 0.16302)

def c1_calcs(c7: float, T: float, beam: float, iE: float) -> float:
    """Calculate the regression-derived coefficient c1, based on c7, the draught and beam, and the half-angle of entrance.

    :param c7: The c7 coefficient.
    :type c7: float
    :param T: Draft of the ship.
    :type T: float
    :param beam: Beam of the ship.
    :type beam: float
    :param iE: Half angle of entrance.
    :type iE: float
    :returns: The c1 coefficient.
    :rtype: float
    """
    return 2223105 * c7 ** 3.78613 * (T / beam) ** 1.07961 * (90 - iE) ** -1.37565

def froude_length_calcs(v: float, length: float, G: float = 9.81) -> float:
    """Calculate the froude number given the velocity, gravity, and waterline length.

    :param v: Ship speed.
    :type v: float
    :param length: Waterline length of the ship.
    :type length: float
    :param G: Acceleration due to gravity, defaults to 9.81.
    :type G: float, optional
    :returns: Froude number.
    :rtype: float
    """
    return v / math.sqrt(length * G)

def lamdba_calcs(length: float, beam: float, cP: float) -> float:
    """Calculate lambda, a regression-derived coefficient relating to wave resistance.

    :param length: Waterline length of the ship.
    :type length: float
    :param beam: Beam of the ship.
    :type beam: float
    :param cP: Prismatic coefficient.
    :type cP: float
    :returns: The lambda coefficient.
    :rtype: float
    """
    if (length / beam) < 12: return 1.446 * cP - 0.03 * length / beam
    else: return 1.446 * cP - 0.36

def c16_calcs(cP: float) -> float:
    """Calculate c16, a regression-derived coefficient based on the prismatic coefficient.

    :param cP: Prismatic coefficient.
    :type cP: float
    :returns: The c16 coefficient.
    :rtype: float
    """
    if cP <= 0.80: return 8.07981 * cP - 13.8763 * cP ** 2 + 6.984388 * cP ** 3
    else: return 1.73014 - 0.7067

def m1_calcs(length: float, T: float, nabla: float, beam: float, c16: float) -> float:
    """Calculate m1, a regression-derived parameter relating to wave resistances.

    :param length: Waterline length of the ship.
    :type length: float
    :param T: Draft of the ship.
    :type T: float
    :param nabla: Moulded displacement volume.
    :type nabla: float
    :param beam: Beam of the ship.
    :type beam: float
    :param c16: The c16 coefficient.
    :type c16: float
    :returns: The m1 parameter.
    :rtype: float
    """
    return 0.0140407 * length / T - 1.75254 * nabla ** (1/3) / length - 4.79323 * beam / length - c16

def c15_calcs(length: float, nabla: float) -> float:
    """Calculate c15, a regression-derived coefficient based on the ratio of the length cubed and nabla.

    :param length: Waterline length of the ship.
    :type length: float
    :param nabla: Moulded displacement volume.
    :type nabla: float
    :returns: The c15 coefficient.
    :rtype: float
    """
    length_nabla = (length ** 3) / nabla
    if length_nabla < 512: return -1.69385
    elif length_nabla > 1727: return 0
    else: return -1.69385 + (length / nabla ** (1.3) - 8.0)/2.36

def m2_calcs(c15: float, cP: float, Fn: float) -> float:
    """Calculate m2, a regression-derived parameter relating to wave resistances.

    :param c15: The c15 coefficient.
    :type c15: float
    :param cP: Prismatic coefficient.
    :type cP: float
    :param Fn: Froude number.
    :type Fn: float
    :returns: The m2 parameter.
    :rtype: float
    """
    return c15 * cP ** 2 * math.exp(-0.1 * Fn ** -2)

def rW_calcs(c1: float, c2: float, c5: float, nabla: float, rho: float, G: float, m1: float, Fn: float, m2: float, lambda_w: float) -> float:
    """Calculate wave resistance based on previously calculated parameters and coefficients.

    :param c1: The c1 coefficient.
    :type c1: float
    :param c2: The c2 coefficient.
    :type c2: float
    :param c5: The c5 coefficient.
    :type c5: float
    :param nabla: Moulded displacement volume.
    :type nabla: float
    :param rho: Density of the fluid.
    :type rho: float
    :param G: Acceleration due to gravity.
    :type G: float
    :param m1: The m1 parameter.
    :type m1: float
    :param Fn: Froude number.
    :type Fn: float
    :param m2: The m2 parameter.
    :type m2: float
    :param lambda_w: The lambda coefficient.
    :type lambda_w: float
    :returns: Wave resistance.
    :rtype: float
    """
    d = -0.9
    return c1 * c2 * c5 * nabla * rho * G * math.exp(m1 * Fn ** d + m2 * math.cos(lambda_w * Fn ** -2))

########## BULBOUS BOW RESISTANCE (R_B) CALCS ###############

def rB_calcs(aBT: float, hB: float, G: float, TF: float, v: float, rho: float) -> float:
    """Calculate the resistance of a bulbous bow given the ship speed and dimensions of the bow.

    :param aBT: Cross-sectional area of the bulbous bow.
    :type aBT: float
    :param hB: Height of the center of the bulbous bow above the keel line.
    :type hB: float
    :param G: Acceleration due to gravity.
    :type G: float
    :param TF: Forward draft of the ship.
    :type TF: float
    :param v: Ship speed.
    :type v: float
    :param rho: Density of the fluid.
    :type rho: float
    :returns: Bulbous bow resistance, 0 if there is no bulbous bow.
    :rtype: float
    """
    #measure of emergence of the bow
    p_B = 0.56 * math.sqrt(aBT) / (TF - 1.5 * hB)

    #Froude Number based on immersion of the bow
    Fni = v / math.sqrt(G * (TF - hB - 0.25 * math.sqrt(aBT)) + 0.15 * v ** 2)

    #0 if no bulbous bow 
    if aBT == 0: return 0
    else: return 0.11 * math.exp(-3 * p_B ** -2) * Fni ** 3 * aBT ** 1.5 * rho * G / (1 + Fni ** 2)



########## TRANSOM STERN RESISTANCE (R_TR) CALCS ###############

def rTR_calcs(aT: float, v: float, G: float, beam: float, cWP: float, rho: float) -> float:
    """Calculate the resistance caused by a transom stern based on the dimensions of the ship and the stern.

    :param aT: Immersed area of the transverse transom.
    :type aT: float
    :param v: Ship speed.
    :type v: float
    :param G: Acceleration due to gravity.
    :type G: float
    :param beam: Beam of the ship.
    :type beam: float
    :param cWP: Waterplane coefficient.
    :type cWP: float
    :param rho: Density of the fluid.
    :type rho: float
    :returns: Transom stern resistance, 0 if there is no transom stern.
    :rtype: float
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

def cA_calcs(TF: float, length: float, cB: float, c2: float) -> float:
    """Calculate the correlation allowance coefficient.

    :param TF: Forward draft of the ship.
    :type TF: float
    :param length: Waterline length of the ship.
    :type length: float
    :param cB: Block coefficient.
    :type cB: float
    :param c2: The c2 coefficient.
    :type c2: float
    :returns: Correlation allowance coefficient (cA).
    :rtype: float
    """
    #c4 factors in the tramsom stern
    if TF / length <= 0.04: c4 = TF / length
    else: c4 = 0.04
    #
    return 0.006 * (length + 100) ** -0.16 - 0.00205 + 0.003 * math.sqrt(length / 7.5) * cB ** 4 * c2 * (0.04 - c4)


def rA_calcs(v: float, cA: float, rho: float, S: float) -> float:
    """Calculate the model-ship correlation resistance.

    :param v: Ship speed.
    :type v: float
    :param cA: Correlation allowance coefficient.
    :type cA: float
    :param rho: Density of the fluid.
    :type rho: float
    :param S: Wetted area of the hull.
    :type S: float
    :returns: Model-ship correlation resistance.
    :rtype: float
    """
    return 0.5 * rho * v **2 * S * cA


########## TOTAL RESISTANCE CALCS ###############


def rTotal_calcs(rF: float, formFactor: float, rAPP: float, rW: float, rB: float, rTR: float, rA: float) -> float:
    """Calculate total resistance force on the ship (in newtons) using previously calculated factors.

    rTotal = cF * (formFactor) + rAPP + rW + rB + rTR + rA

    :param rF: Frictional resistance.
    :type rF: float
    :param formFactor: Form factor (1 + k1).
    :type formFactor: float
    :param rAPP: Appendage resistance.
    :type rAPP: float
    :param rW: Wave resistance.
    :type rW: float
    :param rB: Bulbous bow resistance.
    :type rB: float
    :param rTR: Transom stern resistance.
    :type rTR: float
    :param rA: Model-ship correlation resistance.
    :type rA: float
    :returns: Total resistance force, in newtons.
    :rtype: float
    """
    #total resistance (newtons)
    return rF * (formFactor) + rAPP + rW + rB + rTR + rA

def externalPower_calcs(rTotal: float, v: float) -> float:
    """Calculate the external/effective power applied to the ship as a whole.

    This is NOT the shaft power, which is the power applied to the shafts. Based on the
    formula that power equals force times velocity.

    :param rTotal: Total resistance force.
    :type rTotal: float
    :param v: Ship speed.
    :type v: float
    :returns: External/effective power.
    :rtype: float
    """
    return rTotal * v


def cV_calcs(formFactor: float, cF: float, cA: float) -> float:
    """Calculate the viscous coefficient.

    Factors in the friction of the water on the ship as well as the influence of hull
    form on viscous pressure drag.

    :param formFactor: Form factor (1 + k1).
    :type formFactor: float
    :param cF: Coefficient of friction.
    :type cF: float
    :param cA: Correlation allowance coefficient.
    :type cA: float
    :returns: Viscous coefficient.
    :rtype: float
    """
    return formFactor * cF + cA




########## PROPELLER-VARIABLE CALCS (w, t, eta_R) ###############
#any calculations that differ depending on the number of propellors for a ship. 
#Holtrop & Mennen only account for 1 or 2 propellers, 
#so anything after that is assumed to have similar efficiency to a twin-screw arrangement
#also assume aft draught (T_A)is the same as the amidships draught (T)

def c8_calcs(beam: float, T: float, S: float, length: float, dProp: float) -> float:
    """Calculate the propeller-variable coefficient c8.

    :param beam: Beam of the ship.
    :type beam: float
    :param T: Draft of the ship.
    :type T: float
    :param S: Wetted area of the hull.
    :type S: float
    :param length: Waterline length of the ship.
    :type length: float
    :param dProp: Propeller diameter.
    :type dProp: float
    :returns: The c8 coefficient.
    :rtype: float
    """
    if beam / T < 5:
        return beam * S / (length * dProp * T)
    else:
        return S * (7 * beam / T - 25) / (length * dProp * (beam / T - 3))



def c9_calcs(c8: float) -> float:
    """Calculate the propeller-variable coefficient c9.

    :param c8: The c8 coefficient.
    :type c8: float
    :returns: The c9 coefficient.
    :rtype: float
    """
    if c8 < 28: return c8
    else:
        return 32 - 16 / (c8 - 24)



def c10_calcs(beam: float, length: float) -> float:
    """Calculate the regression-derived propeller-variable coefficient c10.

    :param beam: Beam of the ship.
    :type beam: float
    :param length: Waterline length of the ship.
    :type length: float
    :returns: The c10 coefficient.
    :rtype: float
    """
    if length / beam > 5.2:
        return beam / length
    else:
        return 0.25 - 0.003328402 * (beam / length - 0.134615385)



def c11_calcs(T: float, dProp: float) -> float:
    """Calculate the regression-derived propeller-variable coefficient c11.

    :param T: Draft of the ship (at the rear of the ship).
    :type T: float
    :param dProp: Propeller diameter.
    :type dProp: float
    :returns: The c11 coefficient.
    :rtype: float
    """
    if T/dProp < 2:
        return T / dProp
    else:
        return 0.0833333 * (T / dProp) ** 3 + 1.33333



def cP1_calcs(cP: float, lcb: float) -> float:
    """Calculate the regression-derived coefficient cP1.

    A prismatic coefficient variable accounting for longitudinal centre of bouyancy.

    :param cP: Prismatic coefficient.
    :type cP: float
    :param lcb: Longitudinal center of bouyancy.
    :type lcb: float
    :returns: The cP1 coefficient.
    :rtype: float
    """
    return 1.45 * cP - 0.315 - 0.0225 * lcb



def c20_calcs(cStern: float) -> float:
    """Calculate the regression-derived coefficient c20 based on the stern coefficient.

    :param cStern: Stern shape coefficient.
    :type cStern: float
    :returns: The c20 coefficient.
    :rtype: float
    """
    return 1.0 + 0.015 * cStern



def c19_calcs(cB: float, cP: float, cM: float) -> float:
    """Calculate the regression-derived coefficient c19 based on the hull form coefficients of the ship.

    :param cB: Block coefficient.
    :type cB: float
    :param cP: Prismatic coefficient.
    :type cP: float
    :param cM: Midship section coefficient.
    :type cM: float
    :returns: The c19 coefficient.
    :rtype: float
    """
    if cP < 0.7:
        return 0.12997 / (0.95 - cB) - 0.11056 / (0.95 - cP)
    else:
        return 0.18567 / (1.3571 - cM) - 0.71276 + 0.38648 * cP



def wt_calcs(c9: float, c20: float, cV: float, length: float, T: float, c11: float, cP1: float, beam: float, c19: float, dProp: float, cStern: float, cB: float, numPropellers: int, c10: float) -> tuple[float, float]:
    """Calculate the wake fraction (w) and thrust deduction coefficients (t).

    :param c9: The c9 coefficient.
    :type c9: float
    :param c20: The c20 coefficient.
    :type c20: float
    :param cV: Viscous coefficient.
    :type cV: float
    :param length: Waterline length of the ship.
    :type length: float
    :param T: Draft of the ship.
    :type T: float
    :param c11: The c11 coefficient.
    :type c11: float
    :param cP1: The cP1 coefficient.
    :type cP1: float
    :param beam: Beam of the ship.
    :type beam: float
    :param c19: The c19 coefficient.
    :type c19: float
    :param dProp: Propeller diameter.
    :type dProp: float
    :param cStern: Stern shape coefficient.
    :type cStern: float
    :param cB: Block coefficient.
    :type cB: float
    :param numPropellers: Number of propellers.
    :type numPropellers: int
    :param c10: The c10 coefficient.
    :type c10: float
    :returns: A tuple of (wake fraction w, thrust deduction coefficient td).
    :rtype: tuple[float, float]
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

def K_calcs(numPropellers: int) -> float:
    """Return K, a constant based on number of propellers.

    Only determined for single or double screw ships; unreliable for more than 2 screws.

    :param numPropellers: Number of propellers.
    :type numPropellers: int
    :returns: The K constant.
    :rtype: float
    :raises Exception: If numPropellers is less than 1.
    """
    if numPropellers == 1: return 0.2
    elif numPropellers > 1: return 0.1
    else: raise Exception("Ship cannot have less than 1 propeller!")

def hShaft_calcs(T: float, propKeelClearance: float, dProp: float) -> float:
    """Return the depth of the shaft centerline.

    Calculated as draught minus prop radius plus prop keel clearance (at the stern).

    :param T: Draft of the ship.
    :type T: float
    :param propKeelClearance: Clearance between the propeller tip and the keel line.
    :type propKeelClearance: float
    :param dProp: Propeller diameter.
    :type dProp: float
    :returns: Depth of the shaft centerline.
    :rtype: float
    """
    return T - (propKeelClearance + dProp / 2)

def pitch_calcs(v: float, w: float, n: float) -> float:
    """Calculate propeller pitch, or how far the propeller travels in one revolution.

    formula: pitch = velocity * wake field coefficient (to account for propeller slip) / rotational frequency

    :param v: Ship speed.
    :type v: float
    :param w: Wake fraction.
    :type w: float
    :param n: Shaft speed (rotations per second).
    :type n: float
    :returns: Propeller pitch.
    :rtype: float
    """
    return v * (1 - w) / n

def propThrust_calcs(rTotal: float, td: float) -> float:
    """Return the thrust generated by the propellers.

    effective thrust = total resistance
    effective thrust = prop thrust * thrust deduction coefficient (1 - t)
    prop thrust = total resistance / thrust deduction coefficient

    :param rTotal: Total resistance force.
    :type rTotal: float
    :param td: Thrust deduction coefficient.
    :type td: float
    :returns: Propeller thrust.
    :rtype: float
    """
    return rTotal / (1 - td)



def bladeAreaRatio_calcs(K: float, numBlades: int, propThrust: float, dProp: float, rho: float, G: float, hShaft: float) -> float:
    """Calculate the ratio of the area of the blades.

    original equation: bladeAreaRatio = K + (1.3 + 0.3 * Z) * thrust / (dProp ** 2 * (p_o + rho * G * h - p_v))
    modified to replace p_o - p_v with 99047 N/m^2, which holds true for seawater at 15 degrees celsius

    :param K: Constant based on number of propellers.
    :type K: float
    :param numBlades: Number of blades per propeller.
    :type numBlades: int
    :param propThrust: Propeller thrust.
    :type propThrust: float
    :param dProp: Propeller diameter.
    :type dProp: float
    :param rho: Density of the fluid.
    :type rho: float
    :param G: Acceleration due to gravity.
    :type G: float
    :param hShaft: Depth of the shaft centerline.
    :type hShaft: float
    :returns: Blade area ratio.
    :rtype: float
    """
    return  K + (1.3 + 0.3 * numBlades) * propThrust / (dProp ** 2 * (99047 + rho * G * hShaft))


def etaR_calcs(numPropellers: int, bladeAreaRatio: float, cP: float, lcb: float, pitch: float, dProp: float) -> float:
    #calculate eta_R here since it is dependent on A_E/A_O for the case of a single screw ship
    if numPropellers == 1:
        return 0.9922 - 0.05908 * bladeAreaRatio + 0.07424 * (cP - 0.0225 * lcb)
    else:
        return 0.9737 + 0.111 * (cP - 0.0225 * lcb) + 0.06325 * pitch / dProp

########## OPEN WATER PROPELLER EFFICIENCY CALCS (eta_o) ###############

def chordLength_calcs(bladeAreaRatio: float, dProp: float, numBlades: int) -> float:
    """Calculate c_0.75, or the chord length of the propeller.

    :param bladeAreaRatio: Blade area ratio.
    :type bladeAreaRatio: float
    :param dProp: Propeller diameter.
    :type dProp: float
    :param numBlades: Number of blades per propeller.
    :type numBlades: int
    :returns: Chord length at 0.75 radius (c_0.75).
    :rtype: float
    """
    return 2.073 * (bladeAreaRatio) * dProp / numBlades

def tc_calcs(numBlades: int, c_075: float, dProp: float) -> float:
    """Return the thickness-chord length ratio of the propeller.

    :param numBlades: Number of blades per propeller.
    :type numBlades: int
    :param c_075: Chord length at 0.75 radius.
    :type c_075: float
    :param dProp: propeller diameter.
    :type dProp: float
    :returns: Thickness-chord length ratio.
    :rtype: float
    """
    return (0.0185 - 0.00125 * numBlades) * dProp / c_075

def deltaCD_calcs(tc_075: float, c_075: float, k_p: float) -> float:
    """Calculate the difference in drag coefficients of the hull profile section.

    :param tc_075: Thickness-chord length ratio.
    :type tc_075: float
    :param c_075: Chord length at 0.75 radius.
    :type c_075: float
    :param k_p: Propeller roughness.
    :type k_p: float
    :returns: Difference in drag coefficients.
    :rtype: float
    """
    return (2 + 4 * tc_075) * (0.003605 - (1.89 + 1.62 * math.log(c_075 / k_p)) ** -2.5)

def ktb_calcs(propThrust: float, rho: float, dProp: float, n: float) -> float:
    """Calculate the thrust coefficient of a model (K_T_B or K_TM).

    Based on the ITTC 1978 specification formula.

    :param propThrust: Propeller thrust.
    :type propThrust: float
    :param rho: Density of the fluid.
    :type rho: float
    :param dProp: Propeller diameter.
    :type dProp: float
    :param n: Shaft speed (rotations per second).
    :type n: float
    :returns: Model thrust coefficient (K_T_B).
    :rtype: float
    """

    #ITTC 1978 formulas
    return propThrust / (rho * dProp ** 4 * n ** 2)




def kt_calcs(K_T_B: float, delta_CD: float, pitch: float, c_075: float, numBlades: int, dProp: float) -> float:
    """Calculate the ship thrust coefficient using the Holtrop-Mennen formulation.

    :param K_T_B: Model thrust coefficient.
    :type K_T_B: float
    :param delta_CD: Difference in drag coefficients.
    :type delta_CD: float
    :param pitch: Propeller pitch.
    :type pitch: float
    :param c_075: Chord length at 0.75 radius.
    :type c_075: float
    :param numBlades: Number of blades per propeller.
    :type numBlades: int
    :param dProp: Propeller diameter.
    :type dProp: float
    :returns: Ship thrust coefficient (K_T).
    :rtype: float
    """
    return K_T_B + delta_CD * 0.3 * (pitch * c_075 * numBlades) / dProp ** 2


def j_calcs(v: float, w: float, td: float, n: float, dProp: float) -> float:
    """Calculate the advance ratio of the propeller J.

    :param v: Ship speed.
    :type v: float
    :param w: Wake fraction.
    :type w: float
    :param td: Thrust deduction coefficient.
    :type td: float
    :param n: Shaft speed (rotations per second).
    :type n: float
    :param dProp: Propeller diameter.
    :type dProp: float
    :returns: Advance ratio (J).
    :rtype: float
    """
    return v * (1 - w * td) / (n * dProp)


#can't figure out torque
#K_Q_B = propTorque / (rho * dProp ** 5 * n ** 2)
#K_Q = K_Q_B - delta_CD * 0.25 * (c_075 * numBlades) / dProp
#give up, can't figure way to derive torque
#eta_o = J * K_T /(2 * math.pi * K_Q)

def cTH_calcs(K_T: float, J: float) -> float:
    """Calculate the thrust coefficient, using the formulas provided by ITTC 1978 proceedings.

    :param K_T: Ship thrust coefficient.
    :type K_T: float
    :param J: Advance ratio.
    :type J: float
    :returns: Thrust coefficient (cTH).
    :rtype: float
    """
    return (K_T / J ** 2) * 8 / math.pi


def etao_calcs(cTH: float, trueEfficiencyCoefficient: float = 0.7) -> float:
    """Calculate eta_o, the propeller efficiency.

    Uses a textbook ideal propeller efficiency formula, which does overestimate propeller
    efficiency. The 0.7 coefficient is derived from fitting to the example ship provided
    by Holtrop & Mennen.

    :param cTH: Thrust coefficient.
    :type cTH: float
    :param trueEfficiencyCoefficient: Correction coefficient applied to the ideal efficiency, defaults to 0.7.
    :type trueEfficiencyCoefficient: float, optional
    :returns: Propeller efficiency (eta_o).
    :rtype: float
    """
    #ideal propeller efficiency / ideal eta_o calculation, is a drastic overestimation of eta_o:
    eta_o = 2 / (1 + math.sqrt(1 + cTH))

    #multiply ideal propeller efficiency by 0.7 to reflect real life conditions
    #still results in optimistic measurements
    #eta_o *= trueEfficiencyCoefficient
    return trueEfficiencyCoefficient * eta_o
########## TOTAL SHAFT POWER CALCS ###############

def shaftPowerCalcs(P_E: float, eta_R: float, eta_o: float, eta_S: float, td: float, w: float) -> float:
    """Calculate the required shaft power given the external power and all the efficiency coefficients.

    :param P_E: External/effective power.
    :type P_E: float
    :param eta_R: Relative rotative efficiency.
    :type eta_R: float
    :param eta_o: Open water propeller efficiency.
    :type eta_o: float
    :param eta_S: Sea efficiency coefficient.
    :type eta_S: float
    :param td: Thrust deduction coefficient.
    :type td: float
    :param w: Wake fraction.
    :type w: float
    :returns: Required shaft power.
    :rtype: float
    """
    return P_E / (eta_R * eta_o * eta_S * (1 - td)/(1 - w))


def HoltropMennenPowerCalculation(length: float, beam: float, T: float, displacementMass: float, v: float,
                                  lcb: float = 0, cM: float = 0.95, sAPP: float = 0, cWP: float = 0.7, aBT: float = 0,
                                  hB: float = 4, aT: float = 0, numPropellers: int = 2, dProp: float = 3.5,
                                  numBlades: int = 3, n: float = 3, propKeelClearance: float = 0.2, trueEfficiencyCoefficient: float = 0.7, flowAppendage = 1.5) -> float:
    """Full calculation of resistance and shaft power. All units are in metric.

    :param length: Waterline length of ship in meters.
    :type length: float
    :param beam: Beam of the ship (at its widest) in meters.
    :type beam: float
    :param T: The average draft/depth of the ship at the given displacement in meters. Program assumes uniform draft.
    :type T: float
    :param displacementMass: The displacement of the ship at a given load in metric tons. Trial measurements are ideal, but other displacements are fine.
    :type displacementMass: float
    :param v: The speed of the ship, in meters per second.
    :type v: float
    :param lcb: Longitudinal center of bouyancy (% of ship's length in front of amidships [0.5 L]), defaults to 0 (perfectly amidships).
    :type lcb: float, optional
    :param cM: Midship section coefficient (midship section area / beam * draft). Merchant ships are 0.9, Bismarck was 0.97, defaults to 0.95.
    :type cM: float, optional
    :param sAPP: Wetted area appendage. Appendages are any underwater structures protruding from the hull, like the rudder and propellers, defaults to 0 if unsure.
    :type sAPP: float, optional
    :param cWP: Water plane coefficient (waterplane area / beam * length). cWP of 1 is a rectangle. Larger ships have a slightly lower cWP than smaller ones. Farraguts were 0.744, post-WWII destroyers were 0.68, Bismarck was 0.66, Minas Geraes was 0.712, defaults to 0.7.
    :type cWP: float, optional
    :param aBT: Cross-sectional area of the bulbous bow in square meters, defaults to 0 for a non-bulbous bow.
    :type aBT: float, optional
    :param hB: Height of the center of the bulbous bow above keel line in meters, defaults to 4. Not used if aBT = 0.
    :type hB: float, optional
    :param aT: Immersed area of the transverse area of the transom at zero speed in square meters, defaults to 0 for no transom stern.
    :type aT: float, optional
    :param numPropellers: Number of shafts driving propellers (thrusters exclusively for steering don't count), defaults to 2.
    :type numPropellers: int, optional
    :param dProp: Average diameter of the propellers in meters, defaults to 3.5.
    :type dProp: float, optional
    :param numBlades: Number of blades on each propeller, defaults to 3.
    :type numBlades: int, optional
    :param n: The average shaft speed in rotations per second (hertz), defaults to 3. HMS Hood was 210 rpm (3.5 s^-1), Bismarck was 265 rpm (4.416 s^-1), USS Massachusetts was 185 rpm (3.083 rps), RMS Carmania was 185 rpm, Minas Geraes was 147 rpm (2.45 rps).
    :type n: float, optional
    :param propKeelClearance: How many meters between the tip of a propeller at its lowest and the keel line, defaults to 0.2.
    :type propKeelClearance: float, optional
    :param trueEfficiencyCoefficient: Correction coefficient applied to the ideal propeller efficiency, defaults to 0.7.
    :type trueEfficiencyCoefficient: float, optional
    :param flowAppendage: weighted average of the appendage resistance factor of all appendages. See (Holthrop and Mennen 167) for how to calculate.
    :type flowAppendage: float, optional
    :returns: Required shaft power (kilowatts).
    :rtype: float
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
    #print(f"Froude Number: {Fn}")
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
    #print(rTotal, rF, formFactor, rAPP, rW, rB, rTR, rA)
    #total resistance deviates 4% from Holthrop + Mennen's example (underestimation of resistance)
    P_E = externalPower_calcs(rTotal, v)
    #print(f"EHP: {round(P_E / 1000)} kW")
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
    tc_075 = tc_calcs(numBlades, c_075, dProp)

    #difference in drag coefficients of the profile section
    delta_CD = deltaCD_calcs(tc_075, c_075, k_p)
    K_T_B = ktb_calcs(propThrust, rho, dProp, n)
    K_T = kt_calcs(K_T_B, delta_CD, pitch, c_075, numBlades, dProp)
    J = j_calcs(v, w, td, n, dProp)
    cTH = cTH_calcs(K_T, J)
    eta_o = etao_calcs(cTH, trueEfficiencyCoefficient)
    shaftPower = shaftPowerCalcs(P_E, eta_R, eta_o, eta_S, td, w)
    return shaftPower

def main() -> None:
    with open("HoltropMennenTest.csv", 'r') as csvfile:
        # creating a csv reader object
        csvreader = csv.reader(csvfile)
        fields = next(csvreader)
        for row in csvreader:
            if row == []:
                break
            for idx in range(1, 18): 
                if row[idx] != '': row[idx] = float(row[idx].replace(',', ''))
            name, length, beam, draft, displacement, speed, numShafts, numBlades = row[0], row[2], row[3], row[4], row[5], row[11], row[13], row[14]
            #speed /= 1.944
            if row[8] != '': cM = row[8]
            else: cM = 0.95
            if row[9] != '': cWP = row[9]
            else: cWP = 0.7
            if row[15] != '': dProp = row[15]
            else: dProp = 3.5
            if row[16] == 1.0: 
                print(f'{name} has bulbous bow!')
                aBT = 10
            else: 
                aBT = 0
            if row[17] != '': propSpeed = row[17]
            else: propSpeed = 3
            shaftPower = HoltropMennenPowerCalculation(length, beam, draft, displacement, speed, cM = cM, cWP = cWP,
                                numPropellers = numShafts, dProp = dProp,
                                numBlades = numBlades, n = propSpeed, aBT = aBT)
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