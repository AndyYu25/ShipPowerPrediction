
from shipSpeedCalc import HoltropMennenPowerCalculation
import csv
from typing import Union,TypeVar
import random
import math
import pandas as pd
def loadCSV(filename: str = "HoltropMennenTest.csv")->tuple[list, list]:
    """
    Loads the csv as a list, after stripping the header row
    Data is formatted as follows from HoltropMennenTest.csv:
    Y: each entry is the actual shaft horsepower of the vessel
    X: a table of a variety of inputs to HoltropMennen Test. The columns are as follows:
    0: length (meters)
    1: beam (meters)
    2: draft (meters)
    3: displacement (mt)
    4: speed (m/s)
    5: # of shafts
    6: # of propeller blades
    7: midship section coefficient
    8: waterplane area coefficient
    9: propeller diameter
    10: propeller rotational speed (hertz)
    """
    with open(filename, 'r') as csvfile:
        reader = csv.reader(csvfile)
        X, Y = [], []
        next(reader) #skip header row
        for row in reader:
            if row == []:
                break
            for idx in range(1, 17): 
                if row[idx] != '': row[idx] = float(row[idx].replace(',', ''))
            Y.append(row[16])
            XRow = [row[2], row[3], row[4], row[5], row[10], row[12], row[13]]
            #speed /= 1.944
            if row[7] != '': XRow.append(row[7])
            else: XRow.append(0.95)
            if row[8] != '': XRow.append(row[8])
            else: XRow.append(0.7)
            if row[14] != '': XRow.append(row[14])
            else: XRow.append(3.5)
            if row[15] != '': XRow.append(row[15])
            else: XRow.append(3)
            X.append(XRow)
    return X, Y

def trainTestSplit(X: list, Y: list, testFraction: int = 0.1)->tuple[list, list, list, list]:
    """
    Given a list splits the list into a training and testing set. 
    testFraction is the fraction of the dataset to allocate to the test set, with a default of 0.1, or 10%, rounded down
    Assumes data is indexed by row.
    returns the training set and the test set in that order as lists
    """
    seed = random.randint(0, 100) #set a uniform seed
    random.Random(seed).shuffle(X) #randomize order for a random partition of the data
    random.Random(seed).shuffle(Y) #randomize output of Y in the same way as X was randomized, as to preserve order

    testSize = int(testFraction * len(X))
    trainX, trainY, testX, testY = X[testSize:], Y[testSize:], X[:testSize], Y[:testSize]
    return trainX, trainY, testX, testY


def tecEvaluation(X: list, Y: list, tecWeights:list)->float:
    """
    Given a list of input parameters X, the expected output Y (shaft horsepower measurements),
    output the normalized root mean squared deviation of the estimate of the power calculations.

    This function serves as the objective function.
    """
    errorSum = 0
    yMax, yMin = max(Y), min(Y)
    for idx in range(len(X)):
        row = X[idx]
        #calculate the trueEfficiencyCoefficient, formulated as a linear equation
        trueEfficiencyCoefficient = tecWeights[0] + tecWeights[1] * row[0] + tecWeights[2] * row[1] + tecWeights[3] * row[2] + tecWeights[4] * row[3] + tecWeights[5] * row[4]
        #calculate power outwith w/ the Holtrop-Mennen Formula
        y_hat = HoltropMennenPowerCalculation(row[0], row[1], row[2], row[3], row[4], numPropellers = row[5], numBlades = row[6],
                                             cM = row[7], cWP = row[8], dProp = row[9], n = row[10], trueEfficiencyCoefficient = trueEfficiencyCoefficient)
        y_hat = y_hat / 1000
        errorSum += (Y[idx] - y_hat) ** 2

    MSE = errorSum / len(X)
    NRMSD = math.sqrt(MSE) / (yMax - yMin)
    return round(NRMSD, 8)


def main():
    X, Y = loadCSV()
    trainX, trainY, testX, testY = trainTestSplit(X, Y, testFraction = 0.1)
    print(tecEvaluation(X, Y, [0.7, 0, 0, 0, 0, 0]))

if __name__ == "__main__":
    main()
