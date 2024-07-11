import math
import csv
import os

class WaginenBSeries:
    def __init__(self):
        """
        Loads the polynomial regression coefficients from saved csv files
        """
        #LOAD THRUST COEFFICIENTS
        with open(os.path.join("WaginenBSeries", "Torque.csv"), 'r') as csvfile:
            csvreader = csv.reader(csvfile)
            self.torqueRegressionCoefficients = [[float(j) for j in i] for i in csvreader]
            #check array dimensions are correct, raise exception if not
            if len(self.torqueRegressionCoefficients) != 47 or len(self.torqueRegressionCoefficients[0]) != 5:
                raise Exception("WaginenBSeries\\torque.csv improperly modified or corrupted! Try reinstalling or manually reentering coefficient values.")
        #LOAD TORQUE COEFFICIENTS    
        with open(os.path.join("WaginenBSeries", "Thrust.csv"), 'r') as csvfile:
            csvreader = csv.reader(csvfile)
            self.thrustRegressionCoefficients = [[float(j) for j in i] for i in csvreader]
            #check array dimensions are correct, raise exception if not
            if len(self.thrustRegressionCoefficients) != 39 or len(self.thrustRegressionCoefficients[0]) != 5:
                raise Exception("WaginenBSeriesThrust.csv improperly modified or corrupted! Try reinstalling or manually reentering coefficient values.")
        #LOAD THRUST CORRECTION COEFFICIENTS
        with open(os.path.join("WaginenBSeries", "ThrustCorrection.csv"), 'r') as csvfile:
            csvreader = csv.reader(csvfile)
            self.thrustCorrectionCoefficients = [[float(j) for j in i] for i in csvreader]
            #check array dimensions are correct, raise exception if not
            if len(self.thrustCorrectionCoefficients) != 9 or len(self.thrustCorrectionCoefficients[0]) != 6:
                raise Exception("WaginenBSeriesThrust.csv improperly modified or corrupted! Try reinstalling or manually reentering coefficient values.")
        #LOAD TORQUE CORRECTION COEFFICIENTS
        with open(os.path.join("WaginenBSeries", "TorqueCorrection.csv"), 'r') as csvfile:
            csvreader = csv.reader(csvfile)
            self.torqueCorrectionCoefficients = [[float(j) for j in i] for i in csvreader]
            #check array dimensions are correct, raise exception if not
            if len(self.torqueCorrectionCoefficients) != 13 or len(self.thrustCorrectionCoefficients[0]) != 6:
                raise Exception("WaginenBSeriesThrust.csv improperly modified or corrupted! Try reinstalling or manually reentering coefficient values.")
    def thrustCoefficient(self, d: float, Z: float, AEA0: float, p: float, J: float)->float:
        """
        Calculates the propeller thrust coefficient K_T of a propeller given basic propeller characteristics.
        These values are normalized for a Reynolds Number of 0.75. 
        Use thrustCoefficientCorrected to get the correct coefficient value for a full-scale reynolds number, or use the ITTC 1978 method.
        d: propeller diameter. Diameter units must be the same as pitch units!
        Z: number of propeller blades
        AEA0: expanded blade area ratio A_E/A_0
        p: propeller pitch. Pitch units must be the same as diameter units!
        J: the advance ratio J
        """
        kT = 0
        for row in self.thrustRegressionCoefficients:
            a_i, b_i, c_i, d_i, e_i = row
            kT += a_i * (J ** b_i) * ((p / d) ** c_i) * (AEA0 ** d_i) * (Z ** e_i)
        return kT
    def torqueCoefficient(self, d: float, Z: int, AEA0: float, p: float, J: float)->float:
        """
        Calculates the propeller torque coefficient K_Q of a propeller given basic propeller characteristics.
        These values are normalized for a Reynolds Number of 0.75. 
        Use torqueCoefficientCorrected to get the correct coefficient value for a full-scale reynolds number, or use the ITTC 1978 method.
        d: propeller diameter. Diameter units must be the same as pitch units!
        Z: number of propeller blades
        p: propeller pitch. Pitch units must be the same as diameter units!
        AEA0: expanded blade area ratio A_E/A_0
        PD: the pitch-diameter ration P/D
        J: the advance ratio J
        """
        kQ = 0
        for row in self.torqueRegressionCoefficients:
            a_i, b_i, c_i, d_i, e_i = row
            kQ += a_i * (J ** b_i) * ((p / d) ** c_i) * (AEA0 ** d_i) * (Z ** e_i)
        return kQ
    def torqueCoefficientCorrected(self, d: float, Z: int, AEA0: float, p: float, J: float, reynolds: float)->float:
        """
        Calculates the non-normalized propeller torque coefficient K_Q of a propeller given basic propeller characteristics and the reynolds number of the ship.
        Use torqueCoefficient to obtain normalized coefficient values to a reynolds number of 2e6.
        d: propeller diameter. Diameter units must be the same as pitch units!
        Z: number of propeller blades
        p: propeller pitch. Pitch units must be the same as diameter units!
        AEA0: expanded blade area ratio A_E/A_0
        PD: the pitch-diameter ration P/D
        J: the advance ratio J
        """        
        deltaKQ = 0
        for row in self.torqueCorrectionCoefficients:
            a_i, b_i, c_i, d_i, e_i, f_i = row
            deltaKQ += a_i * (J ** b_i) * ((p/d) ** c_i) * (AEA0 ** d_i) * (Z ** e_i) * (math.log10(reynolds - 0.301) ** f_i)
        kQ = self.torqueCoefficient(d, Z, AEA0, p, J) 
        return kQ + deltaKQ
    def thrustCoefficientCorrected(self, d: float, Z: int, AEA0: float, p: float, J: float, reynolds: float)->float:
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
            deltaKT += a_i * (J ** b_i) * ((p/d) ** c_i) * (AEA0 ** d_i) * (Z ** e_i) * (math.log10(reynolds - 0.301) ** f_i)
        kT = self.torqueCoefficient(d, Z, AEA0, p, J) 
        return kT + deltaKT

bSeries = WaginenBSeries()
p = 0.7
d = 1
print(bSeries.thrustCoefficient(d, 4, 0.55, p, 0))
print(bSeries.torqueCoefficient(d, 4, 0.55, p, 0))
