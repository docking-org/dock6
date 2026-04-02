#!/usr/bin/env python

import sys
import os

input1 = sys.argv[1]
input2 = sys.argv[2]


if len(sys.argv) != 3:
       print("Incorrect number of command line arguments. This script takes exactly two input files where the first input is the full fragment library containing 10838 entries, the second input is the user created torsion table to be appended into the generic table.")
       sys.exit()

filenames = [input1, input2]

if os.stat(input1).st_size == 0:
        print("File input 1 empty, please review input file specification.")
        sys.exit()

if os.stat(input2).st_size == 0:
        print("File input 2 empty, please review input file specification.")
        sys.exit()

with open('./full_fraglib.dat', 'w+') as outfile:
    for fname in filenames:
        with open(fname) as infile:
            for line in infile:
                outfile.write(line)


with open('./full_sorted_fraglib.dat', 'w+') as outfile2:
   file2 = open('./full_fraglib.dat', 'r+')
   lines = file2.readlines()
   lines.sort()
   for liness in lines:
       outfile2.write(liness)
   file2.close()
   outfile2.close()

lines_seen = set() # holds lines already seen
outfile3 = open('./unique_full_sorted_fraglib.dat', "w+")
for linesss in open('./full_sorted_fraglib.dat', "r"):
    if linesss not in lines_seen: # not a duplicate
        outfile3.write(linesss)
        lines_seen.add(linesss)
outfile3.close()


