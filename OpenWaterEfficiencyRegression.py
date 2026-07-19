
from shipSpeedCalc import HoltropMennenPowerCalculation
import csv
from typing import Union,TypeVar
import random
import math
import pandas as pd
from sklearn.model_selection import train_test_split

def loadCSV(filename: str = "HoltropMennenTest.csv", testFraction: int = 0.1, random_state: int = 42)->tuple[pd.DataFrame, pd.DataFrame]:
    """
    Loads the csv as a dataframe, then split it into two  
    """
    df = pd.read_csv(filename, dtype = {
        'id': str, 
        'year': int,
        'length': float,
        'beam': float,
        'draft': float,
        'displacement': float,
        'cM': float,
        'cWP': float,
        'speed': float,
        'num_shafts': int, 
        'num_blades': float,
        'dprop': float,
        'bulbous_bow': int,
        'prop_speed': float,
        'shp': float
    }, thousands= ",")
    if testFraction == 0:
        return df, None
    train, test = train_test_split(df, test_size = testFraction, random_state=random_state)
    return train, test

def shaftPowerError(X: pd.DataFrame, cM_default, cWP_default, dProp_default, n_default, 
                    aBT_default, sApp, flowAppendage, lcb, hb, propKeelClearance,trueEfficiencyCoefficient)->float:
    """
    Given a list of input parameters X, the expected output Y (shaft horsepower measurements),
    output the normalized root mean squared deviation of the estimate of the power calculations.

    This function serves as the objective function.
    """
    errorSum = 0
    #fill NaN values with default params
    X['cM'] = X['cM'].fillna(cM_default)
    X['cWP'] = X['cWP'].fillna(cWP_default)
    X['dprop'] = X['dprop'].fillna(dProp_default)
    X['prop_speed'] = X['prop_speed'].fillna(n_default)

    #if a ship has a bulbous bow, set bulbous_bow to the aBT hyperparam
    X['bulbous_bow'] = X['bulbous_bow'] * aBT_default 

    for row in X.itertuples(index=True):

        #calculate power outwith w/ the Holtrop-Mennen Formula
        #ignore transom stern for now
        y_hat = HoltropMennenPowerCalculation(
            row.length, 
            row.beam, 
            row.draft, 
            row.displacement, 
            row.speed, 
            lcb = lcb, 
            cM = row.cM, 
            sAPP = sApp, 
            cWP = row.cWP, 
            aBT = row.bulbous_bow, 
            hB = hb, 
            aT = 0, 
            numPropellers = row.num_shafts,
            numBlades = row.num_blades,
            n = row.prop_speed, 
            propKeelClearance=propKeelClearance,
            trueEfficiencyCoefficient = trueEfficiencyCoefficient,
            flowAppendage=flowAppendage)
        y_hat = y_hat / 1000 #convert from watts to kilowatts
        errorSum += abs(row.shp - y_hat) / row.shp

    mean_error = errorSum / len(X)
    #NRMSD = math.sqrt(MSE) / (X['shp'].max() - X['shp'].min())
    return round(mean_error, 8)


def main():
    train, test = loadCSV(testFraction=0)
    print(shaftPowerError(train, 0.95, 0.7, 3.5, 3, 0, 0, 1.5, 0, 4, 0.2, 0.7))
    #trainX, trainY, testX, testY = trainTestSplit(X, Y, testFraction = 0.1)
    #print(tecEvaluation(X, Y, [0.7, 0, 0, 0, 0, 0]))

if __name__ == "__main__":
    main()
