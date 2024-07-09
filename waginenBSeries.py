import math
import csv

class WaginenBSeries:
    def __init__(self):
        """
        Loads the polynomial regression coefficients from saved csv files
        """
        with open("WaginenBSeriesTorque.csv", 'r') as csvfile:
            csvreader = csv.reader(csvfile)
            self.torqueRegressionCoefficients = [i for i in csvreader]
            #check array dimensions are correct, raise exception if not
            if len(self.torqueRegressionCoefficients) != 46 or len(self.torqueRegressionCoefficients[0]) != 5:
                raise Exception("WaginenBSeriesTorque.csv improperly modified or corrupted! Try reinstalling or manually reentering coefficient values.")
            print(self.torqueRegressionCoefficients)
        with open("WaginenBSeriesThrust.csv", 'r') as csvfile:
            csvreader = csv.reader(csvfile)
            self.thrustRegressionCoefficients = [i for i in csvreader]
            if len(self.thrustRegressionCoefficients) != 38 or len(self.thrustRegressionCoefficients[0]) != 5:
                raise Exception("WaginenBSeriesThrust.csv improperly modified or corrupted! Try reinstalling or manually reentering coefficient values.")
            print(self.thrustRegressionCoefficients)
    def thrustCoefficient(self, D: float, Z: float, AEA0: float, PD: float, J: float)->float:
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
        for row in self.torqueRegressionCoefficients:
            a_i, b_i, c_i, d_i, e_i = row
            kT += a_i * (J ** b_i) * (PD ** c_i) * (AEA0 ** d_i) * (Z ** e_i)
        return kT
    def torqueCoefficient(self, d: float, Z: float, AEA0: float, p: float, J: float)->float:
        """
        Calculates the propeller torque coefficient K_Q of a propeller given basic propeller characteristics.
        These values are normalized for a Reynolds Number of 0.75. 
        Use torque CoefficientCorrected to get the correct coefficient value for a full-scale reynolds number, or use the ITTC 1978 method.
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
            kQ += a_i * (J ** b_i) * (p/d ** c_i) * (AEA0 ** d_i) * (Z ** e_i)
        return kQ



