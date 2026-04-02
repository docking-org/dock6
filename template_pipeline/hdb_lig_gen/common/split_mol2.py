#!/usr/bin/env python

#Ryan G. Coleman

import string
import sys

def split_mol2(inputFileName, outputPrefix):
  '''splits mol2 files with lots of molecules into single files'''
  outCount = 0
  outFile = None
  startNew = True
  dockOutput = False  # keeps from splitting twice
  try:
    inputFile = open(inputFileName, 'r')
    for line in inputFile:
      strippedLine = string.strip(line)
      if strippedLine.startswith("##########                 Name:"):
        dockOutput = True  # recognizes this is a mol2 from DOCK, don't split 2x
      if not startNew and \
          strippedLine.startswith("##########                 Name:"):
        startNew = True
        outFile.close()
      if not startNew and not dockOutput and \
          strippedLine.startswith("@<TRIPOS>MOLECULE"):
        startNew = True
        outFile.close()
      if startNew:
        outCount += 1
        try:
          name = string.split(strippedLine)[2]
        except IndexError:  # no name
          name = ''
        outFileName = outputPrefix + "." + string.zfill(outCount, 5) + \
            "." + name + ".mol2"
        outFile = open(outFileName, 'w')
        startNew = False
      outFile.write(line)
  except StopIteration:
    pass  # EOF is fine
  finally:
    inputFile.close()
    try:
      outFile.close()
    except AttributeError:
      pass

if -1 != string.find(sys.argv[0], "split_mol2.py"):
  inputFile = sys.argv[1]
  outputPrefix = sys.argv[2]
  split_mol2(inputFile, outputPrefix)
