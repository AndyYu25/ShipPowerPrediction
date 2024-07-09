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
