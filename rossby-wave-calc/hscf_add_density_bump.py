__author__ = 'jacob'
import csv
import os
import numpy as np
from scipy.integrate import quad
import math
from fortranformat import FortranRecordWriter, FortranRecordReader

r_nought = 100
percent_change = 1.1
jmax2 = 258
jmax1 = 257

with open("fort.93", "r") as inputFile:
    header_line = inputFile.readline()
    entry_count = 0
    found = False
    temp_array = []
    with open("fort.2", "w") as header_file:
        header_file.write(header_line)
    for line in inputFile:
        fortran_reader = FortranRecordReader('8(1PE10.3,2X)')
        temp = fortran_reader.read(line)
        '''
        if entry_count + len(temp) > r_nought and not found:
            value_to_change = temp[r_nought - entry_count]
            value_to_change += percent_change * value_to_change
            temp[r_nought - entry_count] = value_to_change
            print("Found!")
            found = True
        else'''
        entry_count += len(temp)
        temp_array.append(temp)
    writing_array = []
    break_count = 0
    for index, element in enumerate(temp_array):
        for j, item in enumerate(element):
            if temp_array[index][j] != None:
                writing_array.append(temp_array[index][j])
            else:
                if break_count == 0:
                    half_point = int(len(writing_array) / 2)
                    print(half_point)
                    print(len(writing_array))
                    writing_array[half_point] = writing_array[half_point] * percent_change
                with open("fort.2", "a") as outfile:
                    fortran_writer = FortranRecordWriter('8(1PE10.3,2X)')
                    output_text = fortran_writer.write(writing_array)
                    outfile.write(output_text)
                    if break_count == 0:
                        outfile.write("\n")
                        break_count = 1
                    writing_array = []