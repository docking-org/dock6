#!/usr/bin/env python
#module for importing/exporting pdb files
#ryan g. coleman ryangc@mail.med.upenn.edu

import string
import math
import sys
import geometry  # distance function
import unionfind2
import collections

#declare this here
radiiDefault = {
    'C': 1.9, 'O': 1.6, 'N': 1.65, 'P': 1.9, 'S': 1.9, 'H': 0.,
    'F': 0., 'I': 0., 'U': 0., 'A': 0., 'B': 0., 'L': 0., '*': 0., 'Z': 0.,
    'D': 0., 'K': 0., 'M': 0.}
    #F should really be about 1.5 but not in fortran so not here
#note that changing these breaks compatibility with trisrf/meshsrf surface
#generation programs and does not actually affect the radii used in those
#processes. in other words don't change them. DON'T DO IT. it won't change
#the radii used AT ALL, it will just break things.
aminoAcid3Codes = [
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"]
aminoAcidCodes = [
    'A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K',
    'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']
occPlace = 0  # the occupancy is first, then the bfactor
bfacPlace = 1  # bfactor
#list of polar hydrogens to keep.
keepPolarH = {
    # H   ALA
    'ALA': ('H  '),
    # H   ARG
    # HE  ARG
    #HH11 ARG
    #HH12 ARG
    #HH21 ARG
    #HH22 ARG
    'ARG': ('H  ', 'HE ', 'HH11', 'HH12', 'HH21', 'HH22'),
    # H   ASN
    #HD21 ASN
    #HD22 ASN
    'ASN': ('H  ', 'HD21', 'HD22'),
    # H   ASP
    'ASP': ('H  '),
    # H   CYS
    # HG  CYS
    'CYS': ('H  ', 'HG '),
    # H   CYX
    'CYX': ('H  '),
    # H   GLN
    #HE21 GLN
    #HE22 GLN
    'GLN': ('H  ', 'HE21', 'HE22'),
    # H   GLU
    'GLU': ('H  '),
    # H   GLY
    'GLY': ('H  '),
    # H   HIS
    # HD1 HIS (change to HID if only)
    # HE2 HIS (change to HIE if only) (if both change to HIP)
    'HIS': ('H  ', 'HD1', 'HE2'),
    'HIP': ('H  ', 'HD1', 'HE2'),
    'HIE': ('H  ', 'HE2'),
    'HID': ('H  ', 'HD1'),
    # H   ILE
    'ILE': ('H  '),
    # H   LEU
    'LEU': ('H  '),
    # H   LYS
    # HZ1 LYS
    # HZ2 LYS
    # HZ3 LYS
    'LYS': ('H  ', 'HZ1', 'HZ2', 'HZ3'),
    # H   MET
    'MET': ('H  '),
    # H   PHE
    'PHE': ('H  '),
    #nothing for proline
    'PRO': (),
    # H   SER
    # HG  SER
    'SER': ('H  ', 'HG '),
    # H   THR
    # HG1 THR
    'THR': ('H  ', 'HG1'),
    # H   TRP
    # HE1 TRP
    'TRP': ('H  ', 'HE1'),
    # H   TYR
    # HH  TYR
    'TYR': ('H  ', 'HH '),
    # H   VAL
    'VAL': ('H  '),
    #add HEM now, maybe not want to special case them all eventually
    'HEM': (),
    'HOH': ('H01', 'H02'),
    'WAT': ('HW1', 'HW2'),
    'ACE': (),
    #add NME (N-methyl amide) as C-terminal capping group (Flavio Ballante)
    'NME': ('HN2')
    }  # end of list

#sometimes you want to use an external radii file (with extreme caution).
def readRadiiFile(radiiFileName):
  '''reads file in "c           1.90" format and returns map from name to radius
  uses column specific format, res name ignored.
  atom__res_radius_
  01234567890123456789'''
  radiiFile = open(radiiFileName, 'r')
  nameRadius = {}  # store in dictionary
  for line in radiiFile:
    try:
      name = string.strip(line[0:5]).upper()
      radius = float(line[10:16])
      nameRadius[name] = radius
    except ValueError:  # ignore where there isn't a float
      pass
  radiiFile.close()
  return nameRadius

def debugCoords(coords, outName=None):
  '''builds fake pdb with the list of coords given, for debugging.'''
  fakeLine = "ATOM     89  N   LEU A  51     " + \
             " 45.172  18.661  -0.866  1.00  0.00           N"
  fakePdb = pdbData()
  for count in range(len(coords)):
    fakePdb.processLine(fakeLine)
    fakePdb.updateNewCoord(len(fakePdb.rawData)-1, coords[count])
  if outName is not None:
    fakePdb.write(outName)
  return fakePdb

#following are for later features, automatic ANISOU removal and automatic
#nonstandard residue HETATM->ATOM replacement
nonStandardResidues = ['PTR', 'TPO', 'SEP', 'ACE', 'KCX', 'MSE', 'CME', 'CSD']
removeLinesStartingWith = ['ANISOU']

def turnListIntoString(resChainList, separator='+'):
  '''turns a list into a string without spaces, stupid helper'''
  if len(resChainList) > 0:
    resChainStr = ""
    for resChain in resChainList:
      resChainStr += separator + resChain
    returnStr = string.join(string.split(resChainStr), "")
    return returnStr  # do this in case chain is " "
  else:
    return separator  # no residues nearby... this is weird

class pdbData(object):
  '''stores data for a pdb file consisting of atoms'''

  def __init__(
      self, filename=False, hetOnly=False, radiiOverride=None,
      ignoreWaters=True):
    '''default constructor takes a pdb file name as input, other ways later'''
    radii = radiiDefault
    if radiiOverride is not None:
      radii = readRadiiFile(radiiOverride)  # override the defaults
    self.__nonZeroRadiiCount = 0
    self.__nonZeroRadiiXYZ = None
    self.__nonZeroRadiiXYZr = None
    self.__factorsByRC = {}
    self.rawData = []
    self.coords = []
    self.radii = []
    self.charges = []
    self.hydroCharges = []
    self.factors = []
    self.atoms = []
    self.resNums = []
    self.resNames = []
    self.altChars = []
    self.chains = []
    self.modelNums = []  # keeps track of NMR models if present
    self.atomToRaw = {}
    self.rawToAtom = {}
    self.ignoreWaters = ignoreWaters
    if filename:
      pdbFile = open(filename, 'r')
      try:
        modelNum = 0
        for line in pdbFile:
          if (string.find(line, "MODEL", 0, 5) != -1):
            try:
              modelNum = int(line.split()[1])  # [0] is MODEL, [1] is the number
            except IndexError:
              modelNum = 0
          if not hetOnly or line.startswith("HETATM"):
            self.processLine(line, modelNum, radii)
      except StopIteration:
        pass  # read the end of the file
      pdbFile.close()

  def processLines(self, lines):
    '''alternate way to initialize if just given a list of lines'''
    modelNum = 0
    for line in lines:
      if (string.find(line, "MODEL", 0, 5) != -1):
        modelNum = int(line.split()[1])  # [0] is MODEL, [1] is the number
      self.processLine(line, modelNum)

  def processLine(self, line, modelNumber=0, radiiToUse=radiiDefault):
    if (
        (string.find(line, "ATOM", 0, 4) != -1) or
        (string.find(line, "HETATM", 0, 6) != -1)):
      name = line[12:16]
      if name[0] == ' ':  # because apparently hetatm entries can start one col
        name = name[1:]  # before atom entries (for the atom name)
      else:
        try:
          numberNotLetter = int(name[0])
          name = name[1:]  # otherwise would have triggered exception
        except ValueError:
          pass  # first character is not a number
      altChar = line[16]
      resName = line[17:20]
      if not (self.ignoreWaters and name == 'HOH'):
        self.rawData.append(line)
        x = float(line[30:38])
        y = float(line[38:46])
        z = float(line[46:54])
        chain = line[21:22]
        self.coords.append((x, y, z))
        self.atoms.append(name)
        try:
          radius = radiiToUse[name[0]]
        except KeyError:
          radius = 0.  # default is to ignore
        #sometimes hydrogens are formatted badly and don't have H in 1st column
        if radius > 0. and name.count('H') > 0:  # any H's are bad
          radius = 0.
        self.radii.append(radius)
        if radius > 0.:
          self.__nonZeroRadiiCount += 1
        self.modelNums.append(modelNumber)
        factorStrings = (line[55:60], line[61:67])
        try:
          factors = (float(factorStrings[0]), float(factorStrings[1]))
        except ValueError:  # in case no factors or occupancies present
          factors = (0., 0.)  # set to zero for now
        self.resNums.append(int(line[22:26]))
        self.chains.append(chain)
        self.altChars.append(altChar)  # alternate sidechain conf letters
        self.resNames.append(resName)
        self.factors.append(factors)
        self.atomToRaw[len(self.atoms) - 1] = len(self.rawData) - 1
        self.rawToAtom[len(self.rawData) - 1] = len(self.atoms) - 1

  def write(self, filename):
    '''opens the file, writes to it, closes it'''
    fileOut = open(filename, 'w')
    self.outputLines(fileOut)
    fileOut.close()

  def replaceHETATMwithATOM(self):
    '''replaces 'HETATM' with 'ATOM  ' in the rawdataline'''
    for rawDataCount, rawDataLine in enumerate(self.rawData):
      if rawDataLine[0:6] == 'HETATM':
        newDataLine = 'ATOM  ' + rawDataLine[6:]
        self.rawData[rawDataCount] = newDataLine

  def outputLines(self, fileOut):
    '''writes the lines to an already open file'''
    for rawDataLine in self.rawData:
      if rawDataLine:
        fileOut.write(rawDataLine)
        if not rawDataLine.endswith('\n'):
          fileOut.write('\n')

  def getOutputLines(self):
    '''same output as outputLines, but to a returned list of strings'''
    returnedList = []
    for rawDataLine in self.rawData:
      if rawDataLine:
        if not rawDataLine.endswith('\n'):
          returnedList.append(rawDataLine + '\n')
        else:
          returnedList.append(rawDataLine)
    return returnedList

  def copy(self):
    newPdb = pdbData(False)
    newPdb.__nonZeroRadiiCount = 0
    newPdb.__nonZeroRadiiXYZ = False
    newPdb.__nonZeroRadiiXYZr = False
    newPdb.__factorsByRC = {}
    newPdb.rawData = self.rawData[:]  # copy everything
    newPdb.coords = self.coords[:]
    newPdb.atoms = self.atoms[:]
    newPdb.radii = self.radii[:]
    newPdb.charges = self.charges[:]
    newPdb.hydroCharges = self.hydroCharges[:]
    newPdb.factors = self.factors[:]
    newPdb.resNums = self.resNums[:]
    newPdb.chains = self.chains[:]
    newPdb.resNames = self.resNames[:]
    newPdb.altChars = self.altChars[:]
    newPdb.modelNums = self.modelNums[:]
    newPdb.atomToRaw = self.atomToRaw.copy()
    newPdb.rawToAtom = self.rawToAtom.copy()
    return newPdb

  def getOrderedRawIndices(self):
    '''gets a list of the ordered raw data indices'''
    indices = list(self.atomToRaw.values())
    indices.sort()
    return indices

  def updateNewXYZ(self, rawDataIndex, newX, newY, newZ):
    self.coords[self.rawToAtom[rawDataIndex]] = (newX, newY, newZ)
    newData = self.rawData[rawDataIndex][:31]
    for data in (newX, newY, newZ):
      tempData = "%+6.3f" % data
      tempInt = float(tempData)
      if math.floor(abs(tempInt)) < 10:
        newData += " "
      if math.floor(abs(tempInt)) < 100:
        newData += "%+6.3f" % data + " "
      elif math.floor(abs(tempInt)) < 1000:  # one less decimal place
        newData += "%+6.2f" % data + " "
      else:  # two less decimal places... if this happens the pdb is huge
        newData += "%+6.1f" % data + " "
    newData += " " + self.rawData[rawDataIndex][56:]
    newDataNoPlus = string.replace(newData, "+", " ")
    self.rawData[rawDataIndex] = newDataNoPlus

  def updateNewCoord(self, rawDataIndex, newCoord):
    self.updateNewXYZ(rawDataIndex, newCoord[0], newCoord[1], newCoord[2])

  def clearFactor(self, rawDataIndex, whichFactor=bfacPlace):
    oldFactors = self.factors[rawDataIndex]
    if whichFactor == bfacPlace:
      newFactors = (oldFactors[occPlace], 0.)
    elif whichFactor == occPlace:
      newFactors = (0., oldFactors[bfacPlace])
    self.updateFactors(rawDataIndex, newFactors)

  def updateFactors(self, rawDataIndex, newFactors):
    '''changes the bfactor and occupancy column. now uses longer format for
    occupancy than might be standard, to get more accurate charges'''
    self.factors[rawDataIndex] = newFactors
    #print self.rawData[rawDataIndex]  # start
    newData = "%+2.2f %+2.4f " % (newFactors[0], newFactors[1])
    newData = self.rawData[rawDataIndex][:55] + newData + "\n"
    newDataNoPlus = string.replace(newData, "+", " ")
    self.rawData[rawDataIndex] = newDataNoPlus

  def setRadiiToFactors(self):
    '''use each atom's radii as the bfactor. leave occupancy alone'''
    for rawDataIndex in range(len(self.rawData)):
      newFactors = (
          self.radii[rawDataIndex], self.factors[rawDataIndex][occPlace])
      self.updateFactors(rawDataIndex, newFactors)
    #nothing to return, factors have been modified in place

  def setChargeToFactors(self, chargeZero=False):
    '''use each atom's charge, set to the occupancy column.'''
    for rawDataIndex in range(len(self.rawData)):
      if not chargeZero:
        newFactors = (
            self.factors[rawDataIndex][occPlace], self.charges[rawDataIndex])
      else:
        newFactors = (self.factors[rawDataIndex][occPlace], 0.0)
      self.updateFactors(rawDataIndex, newFactors)
    #nothing to return, factors have been modified in place

  def getModelNumbers(self):
    '''gets a list of unique model numbers (from NMR data)'''
    uniqueModNums = []
    for modelNum in self.modelNums:
      if modelNum not in uniqueModNums:
        uniqueModNums.append(modelNum)
    uniqueModNums.sort()
    return uniqueModNums

  def getOneModel(self, modelNumber):
    '''copies the pdb, removes all but one NMR model, returns new pdb'''
    newPdb = self.copy()
    markedForRemoval = []
    for index, thisModNum in enumerate(newPdb.modelNums):
      if thisModNum != modelNumber:
        markedForRemoval.append(index)
    for index in markedForRemoval:
      newPdb.removeLine(newPdb.atomToRaw[index])
    return newPdb

  def getEachModel(self):
    '''finds every model in this pdb file, returns each as a separate pdb
    object'''
    retPdb = []
    for modelNumber in self.getModelNumbers():
      retPdb.append(self.getOneModel(modelNumber))
    return retPdb

  def getEachModelJustXyz(self):
    '''finds every model in this pdb file, returns each as a separate pdb
    xyz list'''
    retPdb = []
    for modelNumber in self.getModelNumbers():
      retPdb.append(self.getOneModel(modelNumber).coords)
    return retPdb

  def getEachModelJustXyz_restrictAtomNames(self, atomName):
    '''finds every model in this pdb file, returns each as a separate pdb
    xyz list, restricts to just the atom names passed in (one atom per res)'''
    retPdb = []
    for modelNumber in self.getModelNumbers():
      #print modelNumber
      modelPdb = self.getOneModel(modelNumber)
      modelPdb.keepAtomMatching(atomName)
      modelLines = modelPdb.getOutputLines()
      newPdb = pdbData(False)
      newPdb.processLines(modelLines)
      retPdb.append(newPdb.coords)
    return retPdb

  def divisiveClusteringPdb(self, maxClusters=30):
    '''gets all the models in this pdb file, does divisive clustering'''
    eachModelXyz = self.getEachModelJustXyz()
    import divisive_clustering
    clusters = divisive_clustering.divisiveClustering(eachModelXyz, maxClusters)
    return clusters

  def divisiveClusteringPdb_justC(self, maxClusters=30):
    '''gets all the models in this pdb file, does divisive clustering.
    since proteins are large, just use the carbon in the backbone to represent
    each residue. other options later maybe.'''
    eachModelXyz = self.getEachModelJustXyz_restrictAtomNames("C  ")
    import divisive_clustering
    clusters = divisive_clustering.divisiveClustering(eachModelXyz, maxClusters)
    return clusters

  def getAllResidues(self):
    '''returns 2 lists, one of numbers and one of names, for all unique ress'''
    nums, names = [], []  # cumulators
    curNum, curName = None, None
    for index, thisResNum in enumerate(self.resNums):
      thisResName = self.resNames[index]
      if curNum is None or curNum != thisResNum or curName != thisResName:
        curNum = thisResNum
        curName = thisResName
        nums.append(curNum)
        names.append(curName)
    return nums, names

  def getAllResidueNumbers(self):
    '''returns a list of all residue numbers'''
    return self.getAllResidues()[0]

  def getAllResidueNames(self):
    '''returns a list of all residue names'''
    return self.getAllResidues()[1]

  def clearFactorsResidues(self, residueNumbers, matching=True):
    '''copy pdb, for each residue in the list, remove the bfactor column.
    if matching is false, then remove the bfactor if the residue is not in list
    return new pdb'''
    newPdb = self.copy()
    markedForRemoval = []
    for index, thisResNum in enumerate(newPdb.resNums):
      if matching and (thisResNum not in residueNumbers):
        markedForRemoval.append(index)
      elif not matching and (thisResNum in residueNumbers):
        markedForRemoval.append(index)
    for index in markedForRemoval:
      newPdb.clearFactor(newPdb.atomToRaw[index])
    return newPdb

  def getListResidues(self, residueNumbers):
    '''copies the pdb, removes all residues not in the list, returns new pdb'''
    newPdb = self.copy()
    markedForRemoval = []
    for index, thisResNum in enumerate(newPdb.resNums):
      if thisResNum not in residueNumbers:
        markedForRemoval.append(index)
    for index in markedForRemoval:
      newPdb.removeLine(newPdb.atomToRaw[index])
    return newPdb

  def getListResiduesChains(self, residueChainNumbers):
    '''copies the pdb, removes all residues not in the list, returns new pdb'''
    resChainCheck = []
    for resChainNum in residueChainNumbers:
      resChainCheck.append(resChainNum + '+')  # suffix of + required for match
    newPdb = self.copy()
    markedForRemoval = []
    for index, thisResNum in enumerate(newPdb.resNums):
      chain = newPdb.chains[index]
      thisCheck = string.strip('+' + str(thisResNum) + str(chain)) + '+'
      if thisCheck not in resChainCheck:
        markedForRemoval.append(index)
    for index in markedForRemoval:
      newPdb.removeLine(newPdb.atomToRaw[index])
    return newPdb

  def getFactorsByResidueChain(self, resNum, chain):
    '''gets all the factors (burial depths) for a given resnum + chain combo'''
    hashString = str(resNum) + str(chain)
    if hashString in self.__factorsByRC:
      if len(self.__factorsByRC[hashString]) > 0:
        return self.__factorsByRC[hashString]
    #otherwise compute it
    factorsRet = []
    for index, thisChain in enumerate(self.chains):
      if thisChain == chain:
        if resNum == self.resNums[index]:
          if self.radii[index] > 0.:
            factorsRet.append(self.factors[index][1])
    self.__factorsByRC[hashString] = factorsRet
    return factorsRet

  def getOneChain(self, chainId):
    '''copies the pdb, removes all but one chain, returns it'''
    newPdb = self.copy()
    markedForRemoval = []
    for index, thisChain in enumerate(newPdb.chains):
      if thisChain != chainId:
        markedForRemoval.append(index)
    for index in markedForRemoval:
      newPdb.removeLine(newPdb.atomToRaw[index])
    return newPdb

  def removeLine(self, rawDataIndex):
    self.rawData[rawDataIndex] = False

  def removeAtomMatching(self, index, character):
    '''for removing all zeta-atoms for instance (call with 1, 'Z')'''
    markedForRemoval = []
    for atomIndex, atom in enumerate(self.atoms):
      if atom[index] == character:
        markedForRemoval.append(self.atomToRaw[atomIndex])
    for rawIndex in markedForRemoval:
      self.removeLine(rawIndex)
    #no return necessary, as self has been modified

  def keepAtomMatching(self, atomName):
    '''for making CA or C one-atom per residue models, call with "C  " or "CA "
    as atomName, only those will be kept, all others discarded.'''
    markedForRemoval = []
    for atomIndex, atom in enumerate(self.atoms):
      if atom != atomName:
        markedForRemoval.append(self.atomToRaw[atomIndex])
    for rawIndex in markedForRemoval:
      self.removeLine(rawIndex)
    #no return necessary, as self has been modified

  def removeAllHydrogens(self, resList=None):
    '''for each residue in the list, remove all the hydrogens. if no list given,
    delete all hydrogens in whole protein.'''
    markedForRemoval = []
    for atomIndex, atom in enumerate(self.atoms):
      if atom[0] == 'H':  # hydrogen atom
        resNum = self.resNums[atomIndex]
        if (resList is None) or (resNum in resList):
          markedForRemoval.append(self.atomToRaw[atomIndex])
    for rawIndex in markedForRemoval:
      self.removeLine(rawIndex)

  def removeApolarHydrogen(self):
    '''for removing all nonpolar hydrogens in a protein. uses keepPolarH as the
    dict of residue->atom names to decide which hydrogens to keep'''
    markedForRemoval = []
    for atomIndex, atom in enumerate(self.atoms):
      if atom[0] == 'H':  # hydrogen atom
        resName = self.resNames[atomIndex]
        if resName not in list(keepPolarH.keys()):
          print(("ERROR: residue name unknown:", resName, self.rawData[atomIndex]))
        else:
          allowedHydrogens = keepPolarH[resName]
          if atom not in allowedHydrogens:
            markedForRemoval.append(self.atomToRaw[atomIndex])
    for rawIndex in markedForRemoval:
      self.removeLine(rawIndex)
    #no return necessary, as self has been modified

  def removeProtonsForCovalentDocking(self, resNum, resName, atomNames):
    '''given a pdb number of the covalent residue, removes its proton(s) to allow for the 
    covalent bond'''
    for atomName in atomNames.split(','):
      atomIndex = self.getIndexByResidueAtom(resNum, resName, atomName)
      self.removeLine(atomIndex)

  def updateOneResidue(self, rawDataIndex, newNumber):
    '''updates the data and raw line for a residue number'''
    self.resNums[self.rawToAtom[rawDataIndex]] = newNumber
    newData = self.rawData[rawDataIndex][:22]
    tempData = "%4d" % newNumber
    newData += tempData + self.rawData[rawDataIndex][26:]
    self.rawData[rawDataIndex] = newData

  def updateOneResidueName(self, rawDataIndex, newName):
    '''updates the data and raw line for a residue name'''
    self.resNames[self.rawToAtom[rawDataIndex]] = newName
    newData = self.rawData[rawDataIndex][:17]
    newData += newName + self.rawData[rawDataIndex][20:]
    self.rawData[rawDataIndex] = newData

  def replaceAltChars(self, newChar):
    '''updates the data and raw line for all alternate characters'''
    for rawDataIndex in range(len(self.rawData)):
      self.altChars[self.rawToAtom[rawDataIndex]] = newChar
      newData = self.rawData[rawDataIndex][:16]
      newData += newChar + self.rawData[rawDataIndex][17:]
      self.rawData[rawDataIndex] = newData

  def deleteInsertionCodes(self):
    '''insertion codes are sometimes added as 61A for a residue num
    instead of just using 62 (because people want the numbering to line up
    with some other numbering). they cause problems for some people, so
    this is for removing them. they are in column 27'''
    for rawDataIndex in range(len(self.rawData)):
      newData = self.rawData[rawDataIndex][:26]
      newData += ' '
      newData += self.rawData[rawDataIndex][27:]
      self.rawData[rawDataIndex] = newData

  def renameHistidines(self):
    '''renames histidines from HIS to HID, HIE, or HIP based on hydrogens'''
    #two passes, first find histidines named HIS, want to check all at once
    residueSets = {}  # from (chain,res number)  to [indices]
    for rawDataIndex in range(len(self.rawData)):
      if self.resNames[rawDataIndex] == "HIS":  # only care about unchanged his
        chainResNum = (self.chains[rawDataIndex], self.resNums[rawDataIndex])
        if chainResNum not in list(residueSets.keys()):
          residueSets[chainResNum] = []
        residueSets[chainResNum].append(rawDataIndex)
    #second pass, decide whether each HIS is HID, HIE, or HIP
    for indexList in list(residueSets.values()):
      hasD, hasE, hasP = False, False, False
      for rawDataIndex in indexList:
        if self.atoms[rawDataIndex] == 'HD1':
          hasD = True
        if self.atoms[rawDataIndex] == 'HE2':
          hasE = True
      if hasD and hasE:  # both protonated
        hasP = True
      for rawDataIndex in indexList:
        name = 'HID'  # default name is HID, HIS is never allowed.
        if hasD:
          name = 'HID'
        if hasE:
          name = 'HIE'
        if hasP:
          name = 'HIP'
        self.updateOneResidueName(rawDataIndex, name)

  def renameCysteines(self):
    '''renames cysteines from CYS to CYX based on hydrogens'''
    # what about CYM (negative)?
    #two passes, first find cysteines named CYS, want to check all at once
    residueSets = {}  # from (chain,res number)  to [indices]
    for rawDataIndex in range(len(self.rawData)):
      if self.resNames[rawDataIndex] == "CYS":  # only care about unchanged his
        chainResNum = (self.chains[rawDataIndex], self.resNums[rawDataIndex])
        if chainResNum not in list(residueSets.keys()):
          residueSets[chainResNum] = []
        residueSets[chainResNum].append(rawDataIndex)
    #second pass, decide whether each CYS is CYS or CYX
    for indexList in list(residueSets.values()):
      cyx = True
      for rawDataIndex in indexList:
        #print self.atoms[rawDataIndex]
        if self.atoms[rawDataIndex] == 'HG ':
          #print '*********HG******'
          cyx = False
      for rawDataIndex in indexList:
        name = 'CYS'  # default name is CYS.
        if cyx:
          name = 'CYX'
        self.updateOneResidueName(rawDataIndex, name)


  def renameAllChains(self, newChain=' '):
    '''renames chains, default to blank'''
    for rawDataIndex in range(len(self.rawData)):
      self.updateOneChain(rawDataIndex, newChain)

  def updateOneChain(self, rawDataIndex, newChain=' '):
    '''updates the data and raw line for a residue number'''
    self.chains[self.rawToAtom[rawDataIndex]] = newChain
    newData = self.rawData[rawDataIndex][:21]
    tempData = newChain
    newData += tempData + self.rawData[rawDataIndex][22:]
    self.rawData[rawDataIndex] = newData

  def renumberResidues(self):
    '''renumbers residues to be unique as the pdb biological unit doesn't
    do this properly amazingly enough'''
    chainsHere = {}
    for chainId in self.chains:
      chainsHere[chainId] = True
    for chainId in list(chainsHere.keys()):
      self.renumberResiduesOneChain(chainId)  # does all work

  def fixChainIds(self):
    '''sometimes people screw up the chain id and make them all the same.
    go through the residues and if they go down in number then start a new
    chain id.'''
    lastResNum = -10000
    lastChain = '-'
    chainsUsed = []
    replacingChains = False
    newChain = False
    for atomIndex, residueNumber in enumerate(self.resNums):
      if residueNumber < lastResNum:
        replacingChains = True
        newChain = chr(max(ord(self.chains[atomIndex]), max(chainsUsed)) + 1)
      if replacingChains:
        self.updateOneChain(atomIndex, newChain)
      if ord(self.chains[atomIndex]) not in chainsUsed:
        chainsUsed.append(ord(self.chains[atomIndex]))
      lastResNum = residueNumber
      lastChain = self.chains[atomIndex]

  def renumberResiduesOneChain(self, chainId):
    '''renumbers residues on one chain'''
    seenResNums, lastRes, highestResNum = {}, -1000, -1000
    for atomIndex, residueNumber in enumerate(self.resNums):
      if self.chains[atomIndex] == chainId:  # otherwise skip
        if residueNumber > highestResNum:
          highestResNum = residueNumber
        if lastRes != residueNumber:
          if residueNumber in seenResNums:  # renumber
            newNum = seenResNums[residueNumber]
            if newNum:
              newNum = highestResNum + 1
              highestResNum += 1
            seenResNums[residueNumber] = newNum
            rawDataIndex = self.atomToRaw[atomIndex]
            self.updateOneResidue(rawDataIndex, newNum)
          else:
            seenResNums[residueNumber] = True  # indicate we've seen it before
          lastRes = residueNumber  # always do this
        elif lastRes == residueNumber:
          if seenResNums[residueNumber]:  # nothing to do actually
            pass
          else:
            newNum = seenResNums[residueNumber]
            rawDataIndex = self.atomToRaw[atomIndex]
            self.updateOneResidue(rawDataIndex, newNum)

  def renumberResiduesAllChains(self, newNums=None):
    '''renumbers residues on one chain'''
    seenResNums, lastRes, highestResNum = {}, -1000, 0
    if newNums is None:
      newNumIndex = None
    else:
      newNumIndex = -1
    for atomIndex, residueNumber in enumerate(self.resNums):
      rawDataIndex = self.atomToRaw[atomIndex]
      if residueNumber != lastRes:
        if newNums is not None:
          newNumIndex += 1
          newNum = newNums[newNumIndex]
        else:
          highestResNum += 1
          newNum = highestResNum
        self.updateOneResidue(rawDataIndex, newNum)
      else:  # same as last time
        if newNums is not None:
          newNum = newNums[newNumIndex]
        else:
          newNum = highestResNum
        self.updateOneResidue(rawDataIndex, newNum)
      lastRes = residueNumber  # always do this

  def getAverageCoords(self, listPdbFileNames):
    '''takes current pdb and adds the data from the list, averages the coords
    and makes a new pdb with the output. absolutely no error checking is done
    to ensure the atoms are in the same order in the files'''
    newPdb = self.copy()
    otherData = []
    for pdbFileName in listPdbFileNames:
      otherData.append(pdbData(pdbFileName))
    for index, coord in enumerate(self.coords):
      badData = 0
      sums = list(coord)  # instead of tuple
      for data in otherData:
        try:
          for sumIndex in range(len(sums)):
            sums[sumIndex] += data.coords[index][sumIndex]
        except IndexError:  # somehow the model doesn't have the same #atoms
          badData += 1
      for sumIndex in range(len(sums)):
        sums[sumIndex] /= (len(otherData) + 1 - badData)  # badData tracked
      newPdb.updateNewXYZ(index, sums[0], sums[1], sums[2])
    return newPdb

  def calcRMSDfile(self, otherFileName, alphas=False):
    '''wrapper for calcRMSD'''
    otherData = pdbData(otherFileName)
    return self.calcRMSD(otherData, alphas=alphas)

  def calcRMSD(self, other, alphas=False):
    '''calculates rmsd of all atoms between self and other, no checking or
    alignment is done'''
    squaredSum = 0.0
    badData = 0  # keeps track of missing atoms (some NMR models wrong)
    for index, coord in enumerate(self.coords):
      if not alphas or self.atoms[index][0:2] == "CA":
        try:
          otherCoord = other.coords[index]
          squaredSum += geometry.distL2Squared(coord, otherCoord)
        except IndexError:
          badData += 1
      else:
        badData += 1  # also keep track of non-carbon alphas when looking for CA
    #print squaredSum, (len(self.coords) - badData)
    squaredSum /= (len(self.coords) - badData)
    return squaredSum**0.5

  def isPointNearAnyHeavyAtom(self, point, tolerance):
    nonZeroRadiiXyz = self.getHeavyAtomXYZ()  # cache if necessary
    toleranceSquared = tolerance
    for coord in nonZeroRadiiXyz:
      if geometry.distL2Squared(coord, point) < toleranceSquared:
        return True
    return False

  def getHeavyAtomXYZ(self):
    '''returns the xyz of the nonzero radius atoms'''
    if self.__nonZeroRadiiXYZ is None:  # recalculate
      self.__nonZeroRadiiXYZ = []
      for index, radius in enumerate(self.radii):
        if radius > 0.:
          self.__nonZeroRadiiXYZ.append(self.coords[index])
    return self.__nonZeroRadiiXYZ

  def getHeavyAtomXYZRadius(self):
    '''returns the x,y,z,radius of the nonzero radius atoms'''
    if self.__nonZeroRadiiXYZr is None:  # recalculate
      self.__nonZeroRadiiXYZr = []
      for index, radius in enumerate(self.radii):
        if radius > 0.:
          newList = list(self.coords[index])
          newList.append(radius)
          self.__nonZeroRadiiXYZr.append(newList)
    return self.__nonZeroRadiiXYZr

  def getHeavyAtomCount(self):
    '''counts and returns number of heavy atoms'''
    if self.__nonZeroRadiiCount == 0:  # recalculate
      for radius in self.radii:
        if radius > 0.:
          self.__nonZeroRadiiCount += 1
    return self.__nonZeroRadiiCount

  def getIndexByResidueAtom(self, resNum, resCode, atomName):
    '''gets the index matching the input data, returns false if no match'''
    atomNameStr = string.strip(atomName)
    for index in range(len(self.atoms)):
      if resNum == self.resNums[index]:
        if resCode == self.resNames[index]:
          if atomNameStr == string.strip(self.atoms[index]):
            return index
    return False  # not found

  def getListResidueNumberChain(self):
    '''outputs a unique list of residue number and chain combos'''
    listR = []
    for index, resNum in enumerate(self.resNums):
      chain = self.chains[index]
      if 0 == len(listR) or (resNum, chain) != listR[-1]:
        listR.append((resNum, chain))
    return listR

  def assignCharges(self, charge):
    '''assigns charges from the charge object, assumes caller really wants em'''
    self.charges = []  # delete old charges if present
    self.hydroCharges = []
    for index, resName in enumerate(self.resNames):
      atomName = self.atoms[index]
      self.charges.append(charge.getCharge(atomName, resName))
      self.hydroCharges.append(charge.getTrinaryCharge(atomName, resName))

  def clusterAtoms(self, distanceCutoff=2.0):
    '''breaks into distinct unions of atoms based on distance cutoff'''
    ligandClusters = unionfind2.unionFind()
    cutoffSquared = distanceCutoff ** 2.  # faster comparisons
    for index, coord in enumerate(self.coords):
      for index2, coord2 in enumerate(self.coords):
        if index2 > index:  # only do comparisons once each
          distBetweenSquared = geometry.distL2Squared(coord, coord2)
          if distBetweenSquared <= cutoffSquared:
            ligandClusters.union(index, index2)
    clusteredLists = ligandClusters.toLists()
    newPdbs = []  # list of pdbData objects to return
    for oneCluster in clusteredLists:
      newPdb = self.copy()
      markedForRemoval = []
      for index in range(len(self.coords)):
        if index not in oneCluster:
          markedForRemoval.append(index)
      for index in markedForRemoval:
        newPdb.removeLine(newPdb.atomToRaw[index])
      newPdbs.append(newPdb)
    return newPdbs

  def getNearbyAtoms(self, pointList, nearbyDistance=0.):
    '''returns the list of line numbers of atoms '''
    lines = {}
    if nearbyDistance > 0.:  # actually do distance cutoff
      nearbyDistanceSquared = nearbyDistance ** 2.
      for pt in pointList:
        tempSet = set()
        for index, coord in enumerate(self.coords):
          distanceBetween = geometry.distL2Squared(pt, coord)
          if distanceBetween < nearbyDistanceSquared:
            tempSet.update([self.atomToRaw[index]])
        tempList = list(tempSet)
        tempList.sort()
        lines[tuple(pt)] = tempList
    else:  # just find closese atom for each pt
      for pt in pointList:
        bestDist, bestPt = 10000000., False
        for index, coord in enumerate(self.coords):
          distanceBetween = geometry.distL2Squared(pt, coord)
          if distanceBetween < bestDist:
            bestPt = self.atomToRaw[index]
            bestDist = distanceBetween
        lines[tuple(pt)] = [bestPt]  # still needs to be a list
    return lines

  def getResidueChainsFromNums(self, atomList):
    '''returns a list of residuenumber+chain strings from the list of atoms'''
    resChainList = []
    for index in atomList:
      residueNumber = self.resNums[index]
      chain = self.chains[index]
      resChain = str(residueNumber) + str(chain)
      if resChain not in resChainList:  # guarantee uniqueness
        resChainList.append(resChain)
    resChainList.sort()
    return resChainList

  def getResidueNamesChainsFromNums(self, atomList):
    '''returns a list residuename+number+chain strings from the list of atoms'''
    resChainList = []
    for index in atomList:
      residueName = self.resNames[index]
      residueNumber = self.resNums[index]
      chain = self.chains[index]
      resChain = str(residueName) + str(residueNumber) + str(chain)
      if resChain not in resChainList:  # guarantee uniqueness
        resChainList.append(resChain)
    resChainList.sort()  # is this useful anymore?
    return resChainList

  def getResidueNamesChains(self):
    '''returns a list residuename+number+chain strings from all atoms'''
    resChainList = []
    for index in range(len(self.resNames)):
      residueName = self.resNames[index]
      residueNumber = self.resNums[index]
      chain = self.chains[index]
      resChain = str(residueName) + str(residueNumber) + str(chain)
      if resChain not in resChainList:  # guarantee uniqueness
        resChainList.append(resChain)
    resChainList.sort()  # is this useful anymore?
    return resChainList

  def getAtomsFromNums(self, atomList):
    '''returns a list of atom number+name+chain string from the list of atoms'''
    atomChainList = []
    for index in atomList:
      atomName = self.atoms[index]
      chain = self.chains[index]
      atomChain = str(atomName) + str(index) + str(chain)
      if atomChain not in atomChainList:  # guarantee uniqueness
        atomChainList.append(atomChain)
    return atomChainList

  def getNearbyResidues(self, pointList, nearbyDistance=0.):
    '''returns a new pdb of residues in this pdb near the points given'''
    nearbyDistanceSquared = nearbyDistance ** 2.
    residuesNearPts = []
    for pt in pointList:
      for index, coord in enumerate(self.coords):
        distanceBetweenSquared = geometry.distL2Squared(pt, coord)
        if distanceBetweenSquared < nearbyDistanceSquared:
          residueNumber = self.resNums[index]
          chain = self.chains[index]
          resChain = str(residueNumber) + str(chain)
          if resChain not in residuesNearPts:  # guarantee uniqueness
            residuesNearPts.append(resChain)
    residuesNearPts.sort()
    return self.getListResiduesChains(residuesNearPts)

  def countChains(self):
    '''counts the chains'''
    chains = []
    for chainId in self.chains:
      if chainId not in chains:
        chains.append(chainId)
    return len(chains)

  def getFasta(self):
    '''returns a string in FASTA format, one letter per residue'''
    fasta = ""
    lastResChain = "0X"
    for index in range(len(self.resNums)):
      residueName = self.resNames[index]
      residueNumber = self.resNums[index]
      chain = self.chains[index]
      resChain = str(residueNumber) + str(chain)
      if resChain != lastResChain:  # new residue
        try:
          position = aminoAcid3Codes.index(residueName)
          fasta += aminoAcidCodes[position]
        except ValueError:
          fasta += "X"  # unknown amino acid
      lastResChain = resChain
    return fasta

  def getOccupancyResidue(self, resnum):
    '''gets the occupancy for one residue number. just use first found.'''
    for atomCount in range(len(self.atoms)):
      if self.resNums[atomCount] == resnum:
        return self.factors[atomCount][occPlace]

  def getOccupancyResidueChain(self, resnum, chainId):
    '''gets the occupancy for one residue number & chain combination.'''
    for atomCount in range(len(self.atoms)):
      if self.resNums[atomCount] == resnum and \
          self.altChars[atomCount] == chainId:
        return self.factors[atomCount][occPlace]

  def isMostOccupiedResidueChain(self, resnum, chainId):
    '''for all resnums, is the chainId provided the most occupied one (True) or
    not, return False then'''
    highestOccupancy, highestCount = 0.0, 0
    for atomCount in range(len(self.atoms)):
      if self.resNums[atomCount] == resnum:
        if self.factors[atomCount][occPlace] > highestOccupancy:
          highestOccupancy, highestCount = self.factors[atomCount][occPlace], \
              atomCount
    if self.altChars[highestCount] == chainId:
      return True
    else:
      return False

  def selectMostOccupied(self, exceptions=[], leaveAlone=[]):
    '''for each residue with alternate positions, pick the most occupied
    position. if equal, break ties starting with the A position.
    exceptions can be a list of residues to pick the least occupied for.
    leaveAlone is a list of residues that aren't touched at all.'''
    sameAtoms = collections.defaultdict(list)  # collect lists of atoms
    for atomCount in range(len(self.atoms)):
      if self.altChars[atomCount] != ' ':  # space means no alternate position
        if self.resNums[atomCount] not in leaveAlone:
          truple = (self.atoms[atomCount], self.resNames[atomCount],
                    self.resNums[atomCount])
          sameAtoms[truple].append(atomCount)
    for atomTruple, atomList in list(sameAtoms.items()):
      normal = True  # means picked most occupied
      if atomTruple[2] in exceptions:
        normal = False  # means picked least occupied for these residues
      mostOccupied = None
      occupancy = 0.0
      if not normal:
        occupancy = 1.0
      for atomCount in atomList:
        if normal and self.factors[atomCount][occPlace] > occupancy:
          occupancy = self.factors[atomCount][occPlace]
          mostOccupied = atomCount
        elif self.factors[atomCount][occPlace] == occupancy:
          if self.altChars[atomCount] < self.altChars[mostOccupied]:  # A<B<C
            occupancy = self.factors[atomCount][occPlace]
            mostOccupied = atomCount
        elif not normal and self.factors[atomCount][occPlace] < occupancy:
          occupancy = self.factors[atomCount][occPlace]
          mostOccupied = atomCount
      #okay, now remove everything but the mostOccupied
      for atomCount in atomList:
        if atomCount != mostOccupied:  # keep this one
          self.removeLine(self.atomToRaw[atomCount])

  def deleteAlternates(self, only=None):
    '''for each sidechain with alternate positions, delete them all.
    if only exists, only delete residues in the list of residue numbers'''
    sameAtoms = collections.defaultdict(list)  # collect lists of atoms
    for atomCount in range(len(self.atoms)):
      if self.altChars[atomCount] != ' ':  # space means no alternate position
        truple = (self.atoms[atomCount], self.resNames[atomCount],
                  self.resNums[atomCount])
        sameAtoms[truple].append(atomCount)
    for atomTruple, atomList in list(sameAtoms.items()):
      if (only is None) or (atomTruple[2] in only):
        #okay, now remove everything
        for atomCount in atomList:
          self.removeLine(self.atomToRaw[atomCount])

  def deleteAllResidues(self, leaveAlone=None):
    '''deletes all the atoms in the protein, except the residues in
    leaveAlone'''
    sameAtoms = collections.defaultdict(list)  # collect lists of atoms
    for atomCount in range(len(self.atoms)):
      if self.resNums[atomCount] not in leaveAlone:
        truple = (self.atoms[atomCount], self.resNames[atomCount],
                  self.resNums[atomCount])
        sameAtoms[truple].append(atomCount)
    for atomTruple, atomList in list(sameAtoms.items()):
      for atomCount in atomList:
        self.removeLine(self.atomToRaw[atomCount])

  def getAltChars(self, residueNumbers=None):
    '''for each residue is residueNumbers, return the list of alt chars
    (alternate conformations) seen.'''
    returnAltChars = set()
    for atomCount in range(len(self.atoms)):
      if self.resNums[atomCount] in residueNumbers:  # only care about these
        if self.altChars[atomCount] != ' ':  # space means no alternate position
          returnAltChars.add(self.altChars[atomCount])
    return list(returnAltChars)

  def selectOneAlt(self, residueNumbers=None, pickAltChar=None):
    '''for each residue is residueNumbers, salvage only the pickAltChar
    alternate conformation, delete other conformations.'''
    sameAtoms = collections.defaultdict(list)  # collect lists of atoms
    for atomCount in range(len(self.atoms)):
      if self.altChars[atomCount] != ' ':  # space means no alternate position
        truple = (self.atoms[atomCount], self.resNames[atomCount],
                  self.resNums[atomCount])
        sameAtoms[truple].append(atomCount)
    for atomTruple, atomList in list(sameAtoms.items()):
      if atomTruple[2] in residueNumbers:  # means actually pick here
        for atomCount in atomList:
          if self.altChars[atomCount] != pickAltChar:
            self.removeLine(self.atomToRaw[atomCount])

  def residueSets(self):
    '''calculates and returns a list of residue sets'''
    residueSets = {}  # from (chain,res number)  to [indices]
    for rawDataIndex in range(len(self.rawData)):
      chainResNum = (self.chains[rawDataIndex], self.resNums[rawDataIndex])
      if chainResNum not in list(residueSets.keys()):
        residueSets[chainResNum] = []
      residueSets[chainResNum].append(rawDataIndex)
    return residueSets

  def alignOntoThis(self, otherPdb):
    '''aligns each residue from othePdb onto each residue in this pdb.
    returns a new pdb with aligned residues. only uses first 4 atoms for
    alignment. i.e. backbone atoms'''
    import kabsch  # only needed here for now, needs numpy so don't want always
    #import kabsch since some people may not have numpy but still want other
    #things to work correctly
    import numpy  # for error handling
    newPdb = otherPdb.copy()
    otherResSets = otherPdb.residueSets()
    selfResSets = self.residueSets()
    for chainResNum, indexList in list(otherResSets.items()):
      #print chainResNum, indexList
      xyzList = []
      for thisEntry in selfResSets[chainResNum]:
        xyzList.append(self.coords[thisEntry])
      otherList = []
      for otherEntry in indexList:
        otherList.append(otherPdb.coords[otherEntry])
      try:
        newCoords = kabsch.kabsch_align_other_notnumpy(
            otherList[:4], xyzList[:4], otherList)
        #print chainResNum, otherList[:4], xyzList[:4], otherList, newCoords
        for newIndex in indexList:
          newCoord = newCoords.pop(0)
          newPdb.updateNewCoord(newIndex, newCoord)
      except numpy.linalg.linalg.LinAlgError:
        print(('error', numpy.linalg.linalg.LinAlgError, chainResNum))
        print("usually this is end/beginning of chain problems")
    return newPdb  # has sidechains now

  def getXYZatomsMatching(self, atomNames=None):
    '''if atomnames is none, return all xyz coords for all atoms.
    if it is a list of strings, match any atom with one of those names
    and return them'''
    xyzs = []
    for atomIndex, atom in enumerate(self.atoms):
      if atomNames is None or string.strip(atom) in atomNames:
        xyzs.append(self.coords[atomIndex])
    return xyzs

  def alignOntoThisWholeProtein(self, otherPdb, atomNames=['CA']):
    '''does a carbon-alpha only all protein alignment. then rotates the entire
    protein (other) onto this one. returns whole rotated protein.
    other atoms can be used by passing None or a list of desired names.'''
    import kabsch  # only needed here for now, needs numpy so don't want always
    #import kabsch since some people may not have numpy but still want other
    #things to work correctly
    import numpy  # for error handling
    newPdb = otherPdb.copy()
    matchedSelf = self.getXYZatomsMatching(atomNames)
    matchedOther = otherPdb.getXYZatomsMatching(atomNames)
    allOther = otherPdb.getXYZatomsMatching(None)
    print((len(matchedSelf), len(matchedOther), len(allOther)))
    try:
      newCoords = kabsch.kabsch_align_other_notnumpy(
          matchedOther, matchedSelf, allOther)
      for newIndex in range(len(newPdb.rawData)):
        newCoord = newCoords.pop(0)
        newPdb.updateNewCoord(newIndex, newCoord)
    except numpy.linalg.linalg.LinAlgError:
      print(('error', numpy.linalg.linalg.LinAlgError, chainResNum))
    return newPdb  # has sidechains now

#reminder: not part of pdb class
def writeLigandPdb(outName, atomXyzList, atomTypeList, radii, atomCharge):
  '''standalone function to write ligands as pdb files. use with qnifft a lot'''
  outFile = open(outName, 'w')
  for atomCount, atomXyz in enumerate(atomXyzList):  # write atoms to pdb file
    atomType = string.split(atomTypeList[atomCount], '.')[0]  # only use this
    line = 'ATOM  00000  ' + atomType
    for oneSpace in range(3 - len(atomType)):
      line += ' '
    line += ' LIG     1     '
    for atomCoord in atomXyz:
      if math.floor(abs(atomCoord)) < 10:
        line += ' '
      line += '%+6.3f ' % atomCoord
    #need to add radius and charge data as well.
    radius = 0.0
    if atomType in radii:
      radius = radii[atomType]
    line += ' %4.2f' % radius  # add the radius
    if atomCharge is not None:
      line += ' %+6.4f' % atomCharge[atomCount]  # add the charge
    else:
      line += ' %+6.4f' % 0.0  # or a zero
    outFile.write(line + '\n')
  outFile.close()

if -1 != string.find(sys.argv[0], "pdb.py"):
  for filename in sys.argv[1:]:
    pdbD = pdbData(filename)
    print((filename, pdbD.getHeavyAtomCount()))