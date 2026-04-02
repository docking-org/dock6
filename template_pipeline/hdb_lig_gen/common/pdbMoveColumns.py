#!/usr/bin/env python

#be_blasti.csh has a horrible bug where metals get extra columns before
# the coordinate data. rather than fix be_blasti, just work around it.

import string
import sys
import os

def pdbMoveColumns(inputPdbName, outputPdbName):
  dots = (34, 42, 50, 57, 63)  # should be the periods in the columns
  deleteFrom = 22  # where to delete spaces from
  outputPdb = open(outputPdbName, 'w')
  inputPdb = open(inputPdbName, 'r')
  for line in inputPdb:
    #if string.find(line, 'ATOM') == 0:
    if ('ATOM' in line.split()[0]):
      problems = 0
      for dot in dots:
        try:
          if line[dot] != '.':  # this is not good
            problems += 1
        except IndexError:  # might not be that long
          pass  # no problem
      if problems == 0:  # great
        outputPdb.write(line)
      elif problems > 1:  # find where they should be
        extra = 1
        while line[dots[0] + extra] != '.' and extra < 20:  # at 20, give up
          extra += 1
        if extra != 20:
          newLine = line[:deleteFrom]
          newLine += line[deleteFrom + extra:]
        outputPdb.write(newLine)
    else:
      outputPdb.write(line)  # write these anyway
  outputPdb.close()
  inputPdb.close()

#if -1 != string.find(sys.argv[0], 'pdbMoveColumns.py'):  # main program
if ("pdbMoveColumns.py" in  sys.argv[0]):
  inputPdbName = sys.argv[1]
  outputPdbName = sys.argv[2]
  pdbMoveColumns(inputPdbName, outputPdbName)
