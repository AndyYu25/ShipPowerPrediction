import math
import csv
import os

class WaginenBSeries:
    def __init__(self):
        """
        Loads the polynomial regression coefficients from saved csv files
        """
        #LOAD THRUST COEFFICIENTS
        with open(os.path.join("WaginenBSeries", "torque.csv"), 'r') as csvfile:
            csvreader = csv.reader(csvfile)
            self.torqueRegressionCoefficients = [[float(j) for j in i] for i in csvreader]
            #check array dimensions are correct, raise exception if not
            if len(self.torqueRegressionCoefficients) != 47 or len(self.torqueRegressionCoefficients[0]) != 5:
                raise Exception("WaginenBSeries\\torque.csv improperly modified or corrupted! Try reinstalling or manually reentering coefficient values.")
        #LOAD TORQUE COEFFICIENTS    
        with open(os.path.join("WaginenBSeries", "thrust.csv"), 'r') as csvfile:
            csvreader = csv.reader(csvfile)
            self.thrustRegressionCoefficients = [[float(j) for j in i] for i in csvreader]
            #check array dimensions are correct, raise exception if not
            if len(self.thrustRegressionCoefficients) != 39 or len(self.thrustRegressionCoefficients[0]) != 5:
                raise Exception("WaginenBSeriesThrust.csv improperly modified or corrupted! Try reinstalling or manually reentering coefficient values.")
        #LOAD THRUST CORRECTION COEFFICIENTS
        with open(os.path.join("WaginenBSeries", "thrustCorrection.csv"), 'r') as csvfile:
            csvreader = csv.reader(csvfile)
            self.thrustCorrectionCoefficients = [[float(j) for j in i] for i in csvreader]
            #check array dimensions are correct, raise exception if not
            if len(self.thrustCorrectionCoefficients) != 9 or len(self.thrustCorrectionCoefficients[0]) != 6:
                raise Exception("WaginenBSeriesThrust.csv improperly modified or corrupted! Try reinstalling or manually reentering coefficient values.")
        #LOAD TORQUE CORRECTION COEFFICIENTS
        with open(os.path.join("WaginenBSeries", "torqueCorrection.csv"), 'r') as csvfile:
            csvreader = csv.reader(csvfile)
            self.torqueCorrectionCoefficients = [[float(j) for j in i] for i in csvreader]
            #check array dimensions are correct, raise exception if not
            if len(self.torqueCorrectionCoefficients) != 13 or len(self.thrustCorrectionCoefficients[0]) != 6:
                raise Exception("WaginenBSeriesThrust.csv improperly modified or corrupted! Try reinstalling or manually reentering coefficient values.")
    def thrustCoefficient(self, pdRatio: float, Z: float, AEA0: float, J: float)->float:
        """
        Calculates the propeller thrust coefficient K_T of a propeller given basic propeller characteristics.
        These values are normalized for a Reynolds Number of 0.75. 
        Use thrustCoefficientCorrected to get the correct coefficient value for a full-scale reynolds number, or use the ITTC 1978 method.
        pdRatio: Pitch-diameter ratio.
        Z: number of propeller blades
        AEA0: expanded blade area ratio A_E/A_0
        p: propeller pitch. Pitch units must be the same as diameter units!
        J: the advance ratio J
        """
        kT = 0
        for row in self.thrustRegressionCoefficients:
            a_i, b_i, c_i, d_i, e_i = row
            kT += a_i * (J ** b_i) * ((pdRatio) ** c_i) * (AEA0 ** d_i) * (Z ** e_i)
        return kT
    def torqueCoefficient(self, pdRatio: float, Z: int, AEA0: float, J: float)->float:
        """
        Calculates the propeller torque coefficient K_Q of a propeller given basic propeller characteristics.
        These values are normalized for a Reynolds Number of 0.75. 
        Use torqueCoefficientCorrected to get the correct coefficient value for a full-scale reynolds number, or use the ITTC 1978 method.
        pdRatio: Pitch-diameter ratio.
        Z: number of propeller blades
        p: propeller pitch. Pitch units must be the same as diameter units!
        AEA0: expanded blade area ratio A_E/A_0
        J: the advance ratio J
        """
        kQ = 0
        for row in self.torqueRegressionCoefficients:
            a_i, b_i, c_i, d_i, e_i = row
            kQ += a_i * (J ** b_i) * ((pdRatio) ** c_i) * (AEA0 ** d_i) * (Z ** e_i)
        return kQ
    def torqueCoefficientCorrected(self, pdRatio: float, Z: int, AEA0: float, J: float, reynolds: float)->float:
        """
        Calculates the non-normalized propeller torque coefficient K_Q of a propeller given basic propeller characteristics and the reynolds number of the ship.
        Use torqueCoefficient to obtain normalized coefficient values to a reynolds number of 2e6.
        pdRatio: propeller diameter. Diameter units must be the same as pitch units!
        Z: number of propeller blades
        AEA0: expanded blade area ratio A_E/A_0
        PD: the pitch-diameter ration P/D
        J: the advance ratio J
        """        
        deltaKQ = 0
        for row in self.torqueCorrectionCoefficients:
            a_i, b_i, c_i, d_i, e_i, f_i = row
            deltaKQ += a_i * (J ** b_i) * ((pdRatio) ** c_i) * (AEA0 ** d_i) * (Z ** e_i) * ((math.log10(reynolds) - 0.301) ** f_i)
        kQ = self.torqueCoefficient(pdRatio, Z, AEA0, J) 
        return kQ + deltaKQ
    def thrustCoefficientCorrected(self, pdRatio: float, Z: int, AEA0: float, J: float, reynolds: float)->float:
        """
        Calculates the non-normalized propeller thrust coefficient K_T of a propeller given basic propeller characteristics and the reynolds number of the ship.
        Use torqueCoefficient to obtain normalized coefficient values to a reynolds number of 2e6.
        d: propeller diameter. Diameter units must be the same as pitch units!
        Z: number of propeller blades
        p: propeller pitch. Pitch units must be the same as diameter units!
        AEA0: expanded blade area ratio A_E/A_0
        PD: the pitch-diameter ration P/D
        J: the advance ratio J
        """        
        deltaKT = 0
        for row in self.thrustCorrectionCoefficients:
            a_i, b_i, c_i, d_i, e_i, f_i = row
            deltaKT += a_i * (J ** b_i) * ((pdRatio) ** c_i) * (AEA0 ** d_i) * (Z ** e_i) * ((math.log10(reynolds) - 0.301) ** f_i)
        kT = self.torqueCoefficient(pdRatio, Z, AEA0, J) 
        return kT + deltaKT

def main():
    """
    Test using ship example provided by Holtrop & Mennen 
    """
    from shipSpeedCalc import j_calcs, reynolds_calcs, deltaCD_calcs, tc_calcs, chordLength_calcs

    bSeries = WaginenBSeries()
    #print(bSeries.thrustRegressionCoefficients)
    v = 12.8611 #velocity, 25 knots in m/s
    w = 0.2584 # wake fraction
    n = 1.6594 #propeller frequency (hertz)
    diameter = 8
    print(p / diameter)
    propBlades = 4
    ae_ao = 0.7393
    j = j_calcs(v, w, n, diameter)
    print("J: ", j)
    KV = 11.8987e-7 #kinematic viscoscity (m^2/s). Default value of 11.8987e-7 for water at 16 celsius
    length = 205
    #reynolds = reynolds_calcs(length, v, KV)
    #reynolds number based on relative velocity at 0.7R
    c_075 = chordLength_calcs(ae_ao, diameter, propBlades)
    
    #wake speed, v_a
    v_a = v * (1 - w)
    v_r = math.sqrt(v_a ** 2 + (0.7 * math.pi * n * diameter) ** 2)
    c_r = diameter * 0.520 * ae_ao  #chord length estimation at 0.7R, 0.520 is X2 for 4 blades
    reynolds = reynolds_calcs(c_r, v_r, KV)
    print(reynolds)
    print("K_T: ", bSeries.thrustCoefficient(pdRatio, propBlades, ae_ao, p, j))
    print("K_Q: ", bSeries.torqueCoefficient(pdRatio, propBlades, ae_ao, p, j))
    print("Corrected thrust & torque coefficients")
    print(bSeries.thrustCoefficientCorrected(pdRatio, propBlades, ae_ao, p, j, reynolds))
    print(bSeries.torqueCoefficientCorrected(pdRatio, propBlades, ae_ao, p, j, reynolds))

    # Holtrop & Mennen correction:
    k_p = 0.00003
    tc = tc_calcs(propBlades, c_075, diameter)
    deltaCD = deltaCD_calcs(tc, c_075, k_p)
    #k_t_ship = bSeries.thrustCoefficientCorrected(diameter, propBlades, ae_ao, p, 0, reynolds) + deltaCD * 0.3 * p * c_075 * propBlades / diameter ** 2
    #k_q_ship = bSeries.torqueCoefficientCorrected(diameter, propBlades, ae_ao, p, 0, reynolds) - deltaCD * 0.25 * c_075 * propBlades / diameter ** 2
    k_t_ship = bSeries.thrustCoefficient(pdRatio, propBlades, ae_ao, p, j)
    k_q_ship = bSeries.torqueCoefficient(pdRatio, propBlades, ae_ao, p, j)
    print("k_t_ship: ", k_t_ship)
    print("k_q_ship: ", k_q_ship)

    eta_o = j * k_t_ship / (2 * math.pi * k_q_ship)
    print(eta_o)

if __name__ == "__main__":
    main()