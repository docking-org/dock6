#!/usr/bin/env python

#Ryan G. Coleman, Brian K. Shoichet Lab
#converts sybyl atom types from the ligand building toolchain to dock atom
#types (beginning of mol2db2). reads from parameter file (default can be written
#by calling program in appropriate way)

import string
import sys

class AtomConverter(object):
  '''converts from sybyl to dock atom types. reads from parameter file if
  present. otherwise use defaults. can also write defaults used to a parameter
  file that can be edited/tweaked by users who desire such a feature.'''

  #these are the default parameters in a dictionary. it maps all sybyl types
  #from page 7 of:
  # http://tripos.com/tripos_resources/fileroot/pdfs/mol2_format.pdf
  #this file is also ~rgc/Documents/mol2_format.pdf
  # to the dock atom types. currently these are the same as in mol2db
  # and correspond to the atom types in the vdw.parms.amb.mindock file
  #also in ~rgc/Documents/vdw.parms.amb.mindock
  #the format is 'sybylname': dock integer type
  #exceptions are for where sybyl has 1 type and dock has 2
  #the one of these presently is 'H': 6 and 'H-C': 7 where the code
  #deduces that carbon-bonded hydrogens get type 7 though it is general enough
  #to do this for any atom type if the user adds lines to the parameter file
  convertTypesDefault = {'C.3': 5,
                         'C.2': 1,
                         'C.ar': 1,
                         'C.1': 1,
                         'N.3': 10,
                         'N.2': 8,
                         'N.1': 8,
                         'O.3': 12,
                         'O.2': 11,
                         'S.3': 14,
                         'N.ar': 8,
                         'P.3': 13,
                         'H': 6,
                         'H-C': 7,
                         'Br': 17,
                         'Cl': 16,
                         'F': 15,
                         'I': 18,
                         'S.2': 14,
                         'N.pl3': 8,
                         'LP': 25,
                         'Na': 19,
                         'K': 19,
                         'Ca': 21,
                         'Li': 20,
                         'Al': 20,
                         'Du': 25,
                         'Du.C': 25,
                         'Si': 24,
                         'N.am': 8,
                         'S.o': 14,
                         'S.O': 14,
                         'S.o2': 14,
                         'S.O2': 14,
                         'N.4': 9,
                         'O.co2': 11,
                         'C.cat': 1,
                         'H.spc': 6,
                         'O.spc': 11,
                         'H.t3p': 6,
                         'O.t3p': 11,
                         'ANY': 25,
                         'HEV': 25,
                         'HET': 25,
                         'HAL': 25,
                         'Mg': 20,
                         'Cr.oh': 25,
                         'Cr.th': 25,
                         'Se': 25,
                         'Fe': 25,
                         'Cu': 25,
                         'Zn': 26,
                         'Sn': 25,
                         'Mo': 25,
                         'Mn': 25,
                         'Co.oh': 25}

  def __init__(self, parameterFileName=None):
    '''reads the parameter file if present, just reconstruct convertTypes'''
    if parameterFileName is not None:
      self.convertTypes = {}  # new dictionary
      parameterFile = open(parameterFileName, 'r')
      try:
        lines = parameterFile.readlines()
      except StopIteration:
        pass   # end of file
      finally:
        for line in lines:
          tokens = line.split()
          self.convertTypes[tokens[0]] = int(tokens[1])
    else:
      self.convertTypes = self.convertTypesDefault
    self.specialKeys = {}
    for aKey in self.convertTypes.keys():
      #if -1 != string.find(aKey, '-'):  # aKey has a - in it
      #  index = string.find(aKey, '-')
      if -1 != aKey.find('-'):  # aKey has a - in it
        index = aKey.find('-')
        justFirstPart = aKey[:index]
        self.specialKeys[aKey] = justFirstPart

  def printParameters(self):
    '''prints to standard out the parameters used in a readable format'''
    atomkeys = self.convertTypes.keys()
    #atomkeys.sort()
    sorted(atomkeys)
    for key in atomkeys:
      print (key, self.convertTypes[key])

  def convertMol2atomNum(self, mol2data, atomNum):
    '''takes a given atom number from the mol2data and returns the dock
    type translated.'''
    actualNum = atomNum - 1
    actualName = mol2data.atomType[actualNum]
    if actualName not in self.specialKeys.values():  # don't worry about bonds
      return self.convertTypes[actualName]
    else:  # have to worry about bonds
      #for key in self.specialKeys.iterkeys():  # check each special bond
      for key in iter(self.specialKeys.keys()):  # check each special bond
        if mol2data.bondedTo(atomNum, key[-1]):  # found a match
          return self.convertTypes[key]  # so return the exception
      return self.convertTypes[actualName]  # nothing special bonded found

#if -1 != string.find(sys.argv[0], "sybyl2dock.py"):
if -1 != sys.argv[0].find("sybyl2dock.py"):
  #if program is called from the command line, assume user wants a copy
  #of the default parameter file written to standard out. this is the
  #only command use. usually this will be imported and run from
  #somewhere else.
  if len(sys.argv) > 1:
    AtomConverter(sys.argv[1]).printParameters()
  else:
    AtomConverter().printParameters()
