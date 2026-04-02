#!/usr/bin/env python2.7

#Ryan G. Coleman
#solvation file output reader. i think these are from amsol.

import string
import sys

class Solv(object):
  '''reads .solv files from amsol output'''

  def __init__(self, solvFileName=None):
    '''reads in the file'''
    self.name = "fake"
    self.charge = []
    self.polarSolv = []
    self.apolarSolv = []
    self.solv = []
    self.surface = []
    if solvFileName is not None:
      solvfile = open(solvFileName, 'r')
      try:
        for line in solvfile:
          tokens = line.split()
          if self.name == "fake":  # first line... different from the others
            self.name = tokens[0]
            self.totalAtoms = int(tokens[1])
            self.totalCharge = float(tokens[2])
            self.totalPolarSolv = float(tokens[3])
            self.totalSurface = float(tokens[4])
            self.totalApolarSolv = float(tokens[5])
            self.totalSolv = float(tokens[6])
          else:  # all the other lines, one per atom
            try:
              self.charge.append(float(tokens[0]))
              self.polarSolv.append(float(tokens[1]))
              self.surface.append(float(tokens[2]))
              self.apolarSolv.append(float(tokens[3]))
              self.solv.append(float(tokens[4]))
            except ValueError:  # same solv file repeated multiple times????
              raise self.MultiSolvException(line)
      except (StopIteration, self.MultiSolvException):
        pass  # eof okay
      solvfile.close()

  class MultiSolvException(Exception):
    '''stupid exception class since sometimes solv files have more than 1 mol'''

    def __init__(self, value):
      self.value = value

    def __str__(self):
      return repr(self.value)

#if -1 != string.find(sys.argv[0], "solv.py"):
if -1 != sys.argv[0].find("solv.py"):
  #if called from commandline, read in all arguments as mol2 files and print
  #out statistics about the molecule or confirm that it was read in who knows
  for solvfile in sys.argv[1:]:
    solvdata = Solv(solvfile)
    print (solvdata.name, sum(solvdata.charge), solvdata.totalCharge)
