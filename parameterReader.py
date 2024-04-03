import sys
import argparse
import re

class ParameterReader:

    def __init__(self):
        self.parser = argparse.ArgumentParser()
        self.parser.add_argument("-d","--data",required=True,dest="data")
        self.parser.add_argument("-s","--split",default=0,required=False,type=float,dest="split")

    def readParameters(self):
        arguments = self.parser.parse_args()
        parameters = {}
        parameters["data"] = arguments.data
        parameters["splitPercentage"] = arguments.split
        parameters["split"] = parameters["splitPercentage"] > 0
        return parameters


