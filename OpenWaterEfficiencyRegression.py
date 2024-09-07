
from shipSpeedCalc import HoltropMennenPowerCalculation
import csv
from typing import Union,TypeVar
import random

def loadCSV(filename: str = "HoltropMennenTest.csv")->list:
    """
    Loads the csv as a list, after stripping the header row
    """
    with open(filename, 'r') as csvfile:
        reader = csv.reader(csvfile)
        next(reader) #strip header row
        dataset = list(reader)
    return dataset

def trainTestSplit(data: list, testFraction: int = 0.1)->tuple[list, list]:
    """
    Given a list splits the list into a training and testing set. 
    testFraction is the fraction of the dataset to allocate to the test set, with a default of 0.1, or 10%
    Assumes data is indexed by row.
    returns the training set and the test set in that order as lists
    """
    random.shuffle(data) #randomize order for a random partition of the data
    testSize = int(testFraction * len(data))

    return data[testSize:], data[:testSize]

def main():
    dataset = loadCSV()
    train, test = trainTestSplit(dataset)
    print(train)
    print(test)

if __name__ == "__main__":
    main()
