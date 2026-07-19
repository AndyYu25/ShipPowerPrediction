
from shipSpeedCalc import HoltropMennenPowerCalculation
import csv
from typing import Union,TypeVar
import random
import math
import pandas as pd
from sklearn.model_selection import train_test_split
import optuna
from functools import partial

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


def objective(trial: optuna.Trial, X: pd.DataFrame):
    """
    wrapper to include optuna parameter suggestions
    """
    cM_default = trial.suggest_float("cM_default", 0.8, 1.1) #midships section coefficient, clamp up to 0.8 to make sure hull is filled out
    cWP_default = trial.suggest_float("cWP_default", 0.5, 1.0) #clamp up to 0.5 for reasonable hull shape
    dProp_default = trial.suggest_float("dProp_default", 2.0, 6.0) #propeller diameter (meters); set Yamato as upper bound for propeller diameter
    n_default = trial.suggest_float("n_default", 2.0, 6.0) #assume propeller spinning is between 120 rpm and 360 rpm
    aBT_default = trial.suggest_float("aBT_default", 0.0, 100.0) #guess aBT is within one order of magnitude of the Holtrop-Mennen example
    sApp = trial.suggest_float("sApp", 0.0, 500.0) #guess sApp is within one OOM of Holtrop-Mennen example
    flowAppendage = trial.suggest_float("flowAppendage", 1.0, 4.0) #appendage area is 1 + k_2, so must be at least 1, and highest k_2 value is 4.0
    lcb = trial.suggest_float("lcb", -0.5, 0.5) #lcb is expressed as percentage of length relative to amidships, so must be between -0.5 and 0.5
    hb = trial.suggest_float("hb", 0.0, 5.0) #hb is height above keel line, so can't exceed the draft
    propKeelClearance = trial.suggest_float("propKeelClearance", 0.0, 5.0) #How many meters between the tip of a propeller at its lowest and the keel line
    trueEfficiencyCoefficient = trial.suggest_float("trueEfficiencyCoefficient", 0.0, 1.0) #is a coefficient between 0 and 1

    return shaftPowerError(X, cM_default, cWP_default, dProp_default, n_default, aBT_default, sApp, flowAppendage, lcb, hb, propKeelClearance,trueEfficiencyCoefficient)

def main():
    train, test = loadCSV(testFraction=0)
    print(shaftPowerError(train, 0.95, 0.7, 3.5, 3, 0, 0, 1.5, 0, 4, 0.2, 0.7))
    #trainX, trainY, testX, testY = trainTestSplit(X, Y, testFraction = 0.1)
    #print(tecEvaluation(X, Y, [0.7, 0, 0, 0, 0, 0]))
    study = optuna.create_study()
    objective_df = partial(objective, X = train)
    num_steps = 200
    study.optimize(objective_df, n_trials = num_steps, show_progress_bar=True)
    print(optuna.importance.get_param_importances(study))

if __name__ == "__main__":
    main()
