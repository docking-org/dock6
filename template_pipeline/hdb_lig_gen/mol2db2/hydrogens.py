#!/usr/bin/env python

#Ryan G. Coleman
#reads in file containing terminal rotatable hydrogen definitions.
#writes default file

import string
import sys
import math
import geometry
import combinatorics
from collections import defaultdict

class Hydrogens(object):
  '''holds parameters that determine what a rotatable terminal hydrogen is
  rules have the format (1, Atom, bond, Atom, bond, Atom, Degrees) which is:
   Atom - Atom name like "C.ar" or "C" which matches all starting with C
   bond - "1", "2", "3" or "ar" or "*" which means any
   Degrees - "120,240" or "180" or "-" (or something else crazy) "-" means don't
             rotate this hydrogen
  OR the format (2, bond, bond, Atom, bond, Atom, Degrees) which is:
   where the two bonds at front are both applied to the same Atom (the first)
  rules are processed in order so you can exclude certain things with -
  all things not specified are - (no rotation)

  '''
  rulesDefault = [(1, "C.ar", "1", "S", "1", "H", "180"),
                  (1, "C.ar", "1", "O", "1", "H", "180"),
                  (1, "C.1", "1", "S", "1", "H", "-"),
                  (1, "C.1", "1", "O", "1", "H", "-"),
                  (1, "C", "1", "S", "1", "H", "120,240"),
                  (1, "C", "1", "O", "1", "H", "120,240"),
                  (2, "2", "2", "N", "1", "H", "-"),
                  (2, "1", "2", "N", "1", "H", "180")]

  def __init__(self, parameterFileName=None):
    '''constructs from defaults or reads from file'''
    if parameterFileName is not None:
      parameterFile = open(parameterFileName, 'r')
      self.rules = []
      try:
        for line in parameterFile:
          #tokens = string.split(line)
          tokens = line.split()
          self.rules.append((
              int(tokens[0]), tokens[1], tokens[2], tokens[3], tokens[4],
              tokens[5], tokens[6]))
      except StopIteration:
        pass  # EOF
    else:  # no parameter file, use defaults
      self.rules = self.rulesDefault

  def printParameters(self):
    '''prints to standard out the parameters used in a readable format'''
    for rule in self.rules:
      strout = ""
      for part in rule:  # don't print out the distance squared term
        strout = strout + str(part)
      print (strout)  # force newline

  def findTerminalHydrogens(self, mol2data):
    '''takes a Mol2 class, finds all atoms that meet the rules. called first.'''
    atomNums = mol2data.atomNum
    #atomNums.sort()
    sorted(atomNums)
    mol2data.hydrogenRotAngles = []
    mol2data.hydrogensToRotate = 0
    mol2data.dihedrals = None  # set in findDihedrals later
    for count, atomNum in enumerate(atomNums):
      atomType = mol2data.atomType[count]
      result = "-"
      for rule in self.rules:
        if 1 == rule[0]:  # type 1 rule (1, "C.1", "1", "O", "1", "H", "-"),
          if -1 != atomType.find(rule[5]):  # -1 means not found
            if mol2data.bondedTo(atomNum, rule[3], 1, rule[4]):
              if mol2data.bondedTo(atomNum, rule[1], 2, rule[2]):
                result = rule[6]
                break  # quit this, don't look at the rest of the rules
        elif 2 == rule[0]:  # type 2 rule
          if -1 != atomType.find(rule[5]):  # -1 means not found
            if mol2data.bondedTo(atomNum, rule[3], 1, rule[4]):
              if mol2data.bondedTo(atomNum, "", 2, rule[2]):
                if mol2data.bondedTo(atomNum, "", 2, rule[2]):
                  result = rule[6]
                  break  # quit this, don't look at the rest of the rules
      #print atomNum, atomType, result
      mol2data.hydrogenRotAngles.append(result)
      if result != "-":
        mol2data.hydrogensToRotate += 1

  def _findDihedrals(self, mol2data):
    '''private function called from both rotate and reset that finds the atom
    numbers to use for dihedral rotations for any non-"-" hydrogen'''
    if mol2data.dihedrals is None:
      mol2data.dihedrals = {}  # maps atom number to dihedral atom numbers
      mol2data.rotAngles = {}
      #these are 4 atoms that are all bonded in series. last is hydrogen
      for count, angles in enumerate(mol2data.hydrogenRotAngles):
        if angles != "-":  # don't care about the ones that can't rotate/reset
          dihedral = [-1, -1, -1, -1]
          hydrogenNum = mol2data.atomNum[count]
          dihedral[3] = hydrogenNum
          dihedral[2] = mol2data.bondedTo(hydrogenNum, "", 1, None, True)[1]
          dihedral[1] = mol2data.bondedTo(hydrogenNum, "", 2, None, True)[1]
          dihedral[0] = mol2data.bondedTo(hydrogenNum, "", 3, None, True)[1]
          mol2data.dihedrals[hydrogenNum] = dihedral
          mol2data.rotAngles[hydrogenNum] = []
          #for tempAngle in string.split(angles, ","):
          for tempAngle in angles.split(","):
            mol2data.rotAngles[hydrogenNum].append(float(tempAngle))

  def _getCurDihedral(self, atomNum, xyzCount, mol2data):
    '''for a given atomNum, get the dihedral in the mol2'''
    curDihedralNums = mol2data.dihedrals[atomNum]
    curXyz = []
    for curDihedralNum in curDihedralNums:
      curXyz.append(mol2data.getXyz(xyzCount, curDihedralNum))
    dihedral1 = geometry.getDihedralUnited(tuple(curXyz))
    return curXyz, dihedral1

  def resetHydrogens(self, mol2data):
    '''runs after findTerminalHydrogens, resets the rotatable hydrogens to 0
    instead of the potentially bad angle they were set at.'''
    self._findDihedrals(mol2data)
    for xyzCount in range(mol2data.xyzCount):  # iterate over conformations
      for atomIndex, atomNum in enumerate(mol2data.atomNum):
        if atomNum in mol2data.dihedrals.keys() and \
            180. in mol2data.rotAngles[atomNum]:  # only planar get reset
          curXyz, dihedral1 = self._getCurDihedral(atomNum, xyzCount, mol2data)
          #want to find a new rotation angle, 0, and reset
          newTheta = 0 - dihedral1  # radians!!
          #print newTheta
          newHxyz = geometry.rotateAboutLine(
              curXyz[1], curXyz[2], curXyz[3], newTheta)
          #print curXyz[3], newHxyz
          mol2data.atomXyz[xyzCount][atomIndex] = newHxyz
          mol2data.inputHydrogens[xyzCount] = 1  # means reset
          #print curDihedralNums, dihedral1, dihedral1*(180./math.pi)
          #following code checks to make sure everything worked as expected
          #curXyz = []
          #for curDihedralNum in curDihedralNums:
          #  curXyz.append(mol2data.getXyz(xyzCount, curDihedralNum))
          #dihedral2 = geometry.getDihedralUnited(tuple(curXyz))
          #print dihedral1, dihedral2, math.degrees(dihedral2)

  def rotateHydrogens(self, mol2data):
    '''runs after findTerminalHydrogens, rotates the rotatable hydrogens to
    the angles specified by the rules. differs in that it copies everything
    and makes all combinations of hydrogen rotations'''
    self._findDihedrals(mol2data)  # first call in case reset was not done
    rotatedMol2data = mol2data.copy()  # copy this one, put all new in rotated
    rotatedMol2data.atomXyz = []  # clear
    rotatedMol2data.inputEnergy = []  # clear
    rotatedMol2data.inputHydrogens = []  # clear
    rotatedMol2data.origXyzCount = mol2data.xyzCount
    for xyzCount in range(mol2data.xyzCount):  # iterate over orig confs
      indexToCoords = []  # list of lists
      originals = []
      for atomIndex, atomNum in enumerate(mol2data.atomNum):
        if atomNum in mol2data.dihedrals.keys():
          thisCoords = []
          curXyz, dihedral1 = self._getCurDihedral(atomNum, xyzCount, mol2data)
          thisCoords.append(curXyz[3])  # original saved here
          originals.append(curXyz[3])  # also save orig here
          for newAngle in mol2data.rotAngles[atomNum]:
            #want to find a new rotation angle and move to that
            newTheta = math.radians(newAngle)  # convert to radians!!
            #print dihedral1, newAngle, newTheta,
            newHxyz = geometry.rotateAboutLine(
                curXyz[1], curXyz[2], curXyz[3], newTheta)
            thisCoords.append(newHxyz)  # save new ones here
          indexToCoords.append(thisCoords)
      #have to make all possible combinations
      combinations = combinatorics.allCombinations(indexToCoords)
      for aCombo in combinations:
        if aCombo != originals:  # don't need to do anything for the original
          #need to copy inputEnergy and atomXyz, then replace atomXyz
          atomXyzCopy = mol2data.atomXyz[xyzCount][:]
          #replace hydrogen positions here
          replaceIndex = 0
          for atomIndex, atomNum in enumerate(mol2data.atomNum):
            if atomNum in mol2data.dihedrals.keys():
              atomXyzCopy[atomIndex] = aCombo[replaceIndex]  # replace coord
              replaceIndex += 1  # move this counter forward
          rotatedMol2data.atomXyz.append(atomXyzCopy)
          rotatedMol2data.inputEnergy.append(mol2data.inputEnergy[xyzCount])
          rotatedMol2data.inputHydrogens.append(2)  # means rotate
        else:  # aCombo is the original set
          atomXyzCopy = mol2data.atomXyz[xyzCount][:]
          rotatedMol2data.atomXyz.append(atomXyzCopy)
          rotatedMol2data.inputEnergy.append(mol2data.inputEnergy[xyzCount])
          rotatedMol2data.inputHydrogens.append(
              mol2data.inputHydrogens[xyzCount])
    rotatedMol2data.xyzCount = len(rotatedMol2data.atomXyz)  # finally set this
    return rotatedMol2data

  def findAngles(self, mol2data):
    '''runs after findTerminalHydrogens, find the current angles of the
    terminal hydrogens'''
    self._findDihedrals(mol2data)  # first call in case reset was not done
    for xyzCount in range(mol2data.xyzCount):  # iterate over orig confs
      for atomIndex, atomNum in enumerate(mol2data.atomNum):
        if atomNum in mol2data.dihedrals.keys():
          curXyz, dihedral1 = self._getCurDihedral(atomNum, xyzCount, mol2data)
          print ("%+5.10f" % dihedral1)


if __name__ == '__main__':
    if len(sys.argv) > 0:
        import mol2
        mol = mol2.Mol2(sys.argv[1])
        if len(sys.argv) > 2:
            hs = Hydrogens(sys.argv[2])
        else:
            hs = Hydrogens()
        hs.findTerminalHydrogens(mol)
        print (mol.hydrogensToRotate)
    else:
        Hydrogens().printParameters()
