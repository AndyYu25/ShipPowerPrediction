import math


####### INPUT VALUES ############
#user inputs values

length = 205 #length at waterline in meters
beam = 32 #ship beam in meters / moulded breadth
displacement = 51760000 # in kg
T = 10 #average moulded draught/draft (meters)

lcb = -0.75 #longitudinal center of bouynacy (% of ship's length in front of amidships [0.5 L]). Default is 0, or perfectly amidships

cM = 0.98 #midship section coefficient (midship section area / beam * draft). Merchant ships are 0.9, Bismarck was 0.97, high speed destroyer should have 0.8. Default to 0.95
cB = 0.5717 # block coefficient. Generally between 0.45 and 0.65. Generally, larger warships have a higher block coefficient while smaller warships have a lower block coefficient. Default of 0.55

sAPP = 50 #wetted area appendage. Appendages are any underwater structures protruding from the hull, like the rudder and propellors. Put 0 if unsure.

#water plane coefficient (waterplane area / beam * length). cWP of 1 is a rectangle. Larger ships have a slightly lower cWP than smaller ones.
#Farraguts were 0.744, destroyers were 0.68, Bismarck was 0.66. Default of 0.7
cWP = 0.75

aBT = 20 #cross-sectional area of bulbous bow (m^2). 0 for non bulbous bow (default) (m^2)
hB = 4 #height of the center of the bulbous bow above keel line (m). default of 4, not used if aBT = 0

aT = 16 # immersed area of the transverse area of the transom at zero speed. (m^2)
#defaults to 0 for no transom stern

v = 12.8611 #ship speed (m/s)


########## CONSTANTS ###############

G = 9.81 #acceleration due to gravity (m/s^2)
rho  = 997 #density of fluid (water in this case) kg / m^3
KV = 11.8987e-7 #kinematic viscoscity (m^2/s). Default value of 1.1092 for water at 16 celsius

########## FRICTIONAL RESISTANCE (R_F) CALCS ###############


cCrossSection = 0.9 #cross-section coefficient. 1 is a rectangle, 0.5 is a triangle. 0.9 default assuming u-shaped hull
cP = cB / cM #prismatic coefficient calculation

reynolds = length * v / KV #10e6 is unit converions from mm^2 to m^2
cF = 0.075 / (math.log10(reynolds) - 2) ** 2 #cF is coefficient of friction


def c12(T, length): #correct
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

cStern = 10 # 10 for U-shaped section with hogner stern, 0 for normal section, -10 for v-shaped sections
c13 = 1 + 0.003 * cStern

lR = length * ( 1 - cP + 0.06 * cP * lcb / (4 * cP - 1) )

formFactor = c13 * (0.93 + c12(T, length) * (beam / lR) ** 0.92497 * (0.95 * cP) ** -0.521448 * (1 - cP + 0.0225 * lcb) ** 0.6906) # (1 + k1) in original paper


#wetted area of the hull
S = length * (2 * T + beam) * math.sqrt(cM) * (0.453 + 0.4425 * cB - 0.2862 * cM - 0.003467 * beam / T + 0.3696 * cWP) + 2.38 * aBT / cB

rF = 0.5 * rho * cF * S * (v ** 2)

########## APPENDANGE RESISTANCE (R_APP) CALCS ###############



flowAppendage = 1.5 #weighted average approximation. See (Holthrop and Mennen 167) for a more accurate calculation

#resistance from appendages. 1 is the water density
rAPP = 0.5 * rho * (v ** 2) * sAPP * flowAppendage * cF

print(rAPP)
########## Wave RESISTANCE (R_W) CALCS ###############

if (beam / length) < 0.11:
    c7 = 0.229577 * (beam / length) ** (1/3)
elif 0.11 <= beam / length and beam / length <= 0.25:
    c7 = beam / length
elif 0.25 < beam / length:
    c7 = 0.5 - 0.0625 * length / beam

#TF is the forward draft. Ships are assumed to be large and slow enough that the difference between the amidships and forward draft is negligible
TF = T
#influence of bulbous bow on wave resistance
c3 = 0.56 * aBT ** 1.5 / (beam * T * (0.31 * math.sqrt(aBT) + TF - hB))

#c2 is wave resistance coefficient due to bulbous bow. defaults to 1 for non-bulbous bows
if aBT == 0: 
    c2 = 1
else:
    c2 = math.exp(-1.89 * math.sqrt(c3))

#c5 is wave resistance coefficient due to transom sterns. defaults to 1 for non-transom sterns
if aT == 0:
    c5 = 1
else:
    c5 = 1 - 0.8 * aT / (beam * T * cM)

#moulded displacement volume (nabla)
nabla = cB  * length * beam * T

#iE is half angle of entrance (angle of the waterline at the bow with reference to the center plane)
#can be manually inputted, but estimated via regression analysis for user simplicity

iE = 1 + 89 * math.exp(-1 * (length / beam) ** 0.80856 * (1 - cWP) ** 0.30484 * (1 - cP - 0.0225 * lcb) ** 0.6367 
    * (lR / beam) ** 0.34574 * (100 * nabla / length ** 3) ** 0.16302)

c1 = 2223105 * c7 ** 3.78613 * (T / beam) ** 1.07961 * (90 - iE) ** -1.37565

#Froude Number, based on waterline length
Fn = v / math.sqrt(length * G)

if (length / beam) < 12:
    lamdba_w = 1.446 * cP - 0.03 * length / beam
else:
    lamdba_w = 1.446 * cP - 0.36


if cP <= 0.80:
    c16 = 8.07981 * cP - 13.8763 * cP ** 2 + 6.984388 * cP ** 3
else:
    c16 = 1.73014 - 0.7067

m1 = 0.0140407 * length / T - 1.75254 * nabla ** (1/3) / length - 4.79323 * beam / length - c16

length_nabla = (length ** 3) / nabla

d = -0.9


if length_nabla < 512:
    c15 = -1.69385
elif length_nabla > 1727:
    c15 = 0
else:
    c15 = -1.69385 + (length / nabla ** (1.3) - 8.0)/2.36

m2 = c15 * cP ** 2 * math.exp(-0.1 * Fn ** -2)

rW = c1 * c2 * c5 * nabla * rho * G * math.exp(m1 * Fn ** d + m2 * math.cos(lamdba_w * Fn ** -2))

########## BULBOUS BOW RESISTANCE (R_B) CALCS ###############

#measure of emergence of the bow
p_B = 0.56 * math.sqrt(aBT) / (TF - 1.5 * hB)

#Froude Number based on immersion of the bow
Fni = v / math.sqrt(G * (TF - hB - 0.25 * math.sqrt(aBT)) + 0.15 * v ** 2)

#0 if no bulbous bow 
if aBT == 0:
    rB = 0
else:
    rB = 0.11 * math.exp(-3 * p_B ** -2) * Fni ** 3 * aBT ** 1.5 * rho * G / (1 + Fni ** 2)

########## TRANSOM STERN RESISTANCE (R_TR) CALCS ###############
#only run calc if transom stern exists
if aT > 0:
    #Froude Number based on transom immersion
    FnT = v / math.sqrt(2 * G * aT / (beam + beam * cWP))

    if FnT < 5:
        c6 = 0.2 * (1 - 0.2 * FnT)
    else:
        c6 = 0

    rTR = 0.5 * rho * v ** 2 * aT * c6
else:
    rTR = 0

########## MODEL-SHIP CORRELATION RESISTANCE (R_A) CALCS ###############

#c4 factors in the tramsom stern
if TF / length <= 0.04:
    c4 = TF / length
else:
    c4 = 0.04

#correlation allowance coefficient
cA = 0.006 * (length + 100) ** -0.16 - 0.00205 + 0.003 * math.sqrt(length / 7.5) * cB ** 4 * c2 * (0.04 - c4)

rA = 0.5 * rho * v **2 * S * cA

#total resistance (newtons)
#faulty formula + assumptions: disregard resistance from appendages
#rTotal = cF * (formFactor) + rAPP + rW + rB + rTR + rA
rTotal = rF * (formFactor) + rAPP + rW + rB + rTR + rA

#total resistance deviates 4% from provided example (underestimation of resistance)
print(rTotal)

requiredPower = rTotal * v

print(requiredPower)

#sources: overall calculation (Holthrop and Mennen): https://repository.tudelft.nl/islandora/object/uuid:ee370fed-4b4f-4a70-af77-e14c3e692fd4/datastream/OBJ/download
#coefficient of friction + reynolds number calculation: https://repository.tudelft.nl/islandora/object/uuid%3A16d77473-7043-4099-a8c6-bf58f555e2e7
#prismatic coefficient (cP): https://www.nautilusshipping.com/form-coefficient-of-ship
#destroyer ship coefficient numbers: https://www.jstor.org/stable/pdf/44893811.pdf

