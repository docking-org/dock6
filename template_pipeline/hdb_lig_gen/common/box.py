#!/usr/bin/env python
#ryan g. coleman, brian shoichet lab. 2012
#reads/writes box files from DOCK toolchain

import sys
import string
import phi
import math
import combinatorics

def writeBox(minXyz, maxXyz, filename):
  '''writes a "box" file for docking, used in grid construction and not actually
  read by DOCK. allows non-cubic grids.'''
  outFile = open(filename, 'w')
  outFile.write(
      'HEADER    CORNERS OF BOX %8.3f%8.3f%8.3f%8.3f%8.3f%8.3f\n' %
      (minXyz[0], minXyz[1], minXyz[2], maxXyz[0], maxXyz[1], maxXyz[2]))
  outFile.write(
      'REMARK    CENTER (X Y Z) %8.3f%8.3f%8.3f\n' %
      ((maxXyz[0] - minXyz[0]) / 2., (maxXyz[1] - minXyz[1]) / 2.,
       (maxXyz[2] - minXyz[2]) / 2.))
  outFile.write(
      'REMARK    DIMENSIONS (X Y Z) %8.3f%8.3f%8.3f\n' %
      (math.fabs(maxXyz[0] - minXyz[0]), math.fabs(maxXyz[1] - minXyz[1]),
       math.fabs(maxXyz[2] - minXyz[2])))
  names = ['DUA', 'DUB', 'DUC', 'DUD', 'DUE', 'DUF', 'DUG', 'DUH']
  corners = combinatorics.allCombinations(
      [[minXyz[0], maxXyz[0]], [minXyz[1], maxXyz[1]], [minXyz[2], maxXyz[2]]])
  for count in range(8):
    outFile.write(
        'ATOM      %d  %s BOX     1    %8.3f%8.3f%8.3f\n' %
        (count + 1, names[count], corners[count][0], corners[count][1],
         corners[count][2]))
  outFile.write('CONECT    1    2    3    5\n')
  outFile.write('CONECT    2    1    4    6\n')
  outFile.write('CONECT    3    1    4    7\n')
  outFile.write('CONECT    4    2    3    8\n')
  outFile.write('CONECT    5    1    6    7\n')
  outFile.write('CONECT    6    2    5    8\n')
  outFile.write('CONECT    7    3    5    8\n')
  outFile.write('CONECT    8    4    6    7\n')
  outFile.close()

def readBox(filename):
  '''reads a 'box' file from docking, used for constructing other grids.
  since things outside this box can't be scored, we don't need to save that data
  for reading into DOCK. ridiculously allows any odd-numbered gridsize.'''
  corners = []
  center = []
  dimensions = []
  try:
    for line in open(filename, 'r'):
      if line.find("CORNERS") > 0:
        #corners = [float(item) for item in string.split(line)[4:]]
        corners = [float(item) for item in line.split()[4:]]
      elif line.find("CENTER") > 0:
        #center = [float(item) for item in string.split(line)[5:]]
        center = [float(item) for item in line.split()[5:]]
      elif line.find("DIMENSIONS") > 0:
        dimensions = [float(item) for item in line.split()[5:]]
  except StopIteration:
    pass
  return corners, center, dimensions
