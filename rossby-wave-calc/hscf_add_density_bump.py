__author__ = 'jacob'
import csv
import os
import numpy as np
from scipy.integrate import quad
import math
from fortranformat import FortranRecordWriter, FortranRecordReader

r_nought = 100
percent_change = 1.05

with open("fort.93", "r") as inputFile:
    header_line = inputFile.readline()
    print(header_line)
    entry_count = 0
    found = False
    for line in inputFile:
        fortran_reader = FortranRecordReader('8(1PE10.3,2X)')
        temp = fortran_reader.read(line)
        if entry_count + len(temp) > r_nought and not found:
            value_to_change = temp[r_nought - entry_count]
            value_to_change += percent_change * value_to_change
            temp[r_nought - entry_count] = value_to_change
            found = True
        with open("fort.2", "a") as outfile:
            fortran_writer = FortranRecordWriter('8(1PE10.3,2X)')
            print(temp)
            output_text = fortran_writer.write(temp)
            outfile.write(output_text)