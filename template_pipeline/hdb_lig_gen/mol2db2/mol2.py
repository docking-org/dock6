#!/usr/bin/env python2.7

#Ryan G. Coleman
#mol2 reader/writer

import string
import sys
import sybyl2dock
import atom_color_table
import collections
import gzip
import operator
import floydwarshall  # all pairs shortest paths routine
import shortestpaths  # from one point to all others
import geometry  # for distance function
import unionfind2
import divisive_clustering

class Mol2(object):
  '''reads mol2 files into a bunch of lists (for speed instead of objects).
  reads multi-mol2 files containing the same molecule but in diff confs,
  only xyzs are read in separately and no checking is done to ensure the
  molecule is the same. (garbage in garbage out).'''

  def blankNew(self):
    '''create and blank everything'''
    self.name = "fake"  # at some point these need to be capped at 9 characters
    self.protName = "fake"  # yeah this too, 9 characters only
    self.atomNum = []
    self.atomName = []  # all this stuff might need to get optimized
    self.atomXyz = []  # kept for every conformation
    self.inputEnergy = []  # kept for every conformation
    self.inputHydrogens = []  # for every conf. 0 means input. 1 means reset.
                              #2 means rotated. 3 means mix-n-match so who knows
    self.atomType = []
    self.atomCharge = [] #read in but overriden by solv data
    self.atomBonds = []
    self.bondNum = []
    self.bondStart = []
    self.bondEnd = []
    self.bondType = []
    self.xyzCount = -1
    self.origXyzCount = -1
    self.smiles = "fake"
    self.longname = "fake"
    self.bondDists = None
    self.rmsdTable = None

  def __init__(self):
    '''makes a totally fake version, for copying into'''
    self.blankNew()

  def __init__(self, mol2fileName=None, nameFileName=None, mol2text=None):
    '''reads in the file'''
    self.blankNew()
    #read from files/text
    if mol2text is not None:
      for line in mol2text:
        self.processLine(line)
      self.xyzCount += 1
    #read the name.txt file which is one line and is made by the toolchain
    if nameFileName is not None:
      try:
        #name.txt is grepped out from the dict file which has a random format
        namefile = open(nameFileName, 'r')
        firstLine = namefile.readlines()[0]
        tokens = firstLine.split()
        maxsplit = 0
        if 9 == len(tokens[2]) and "P" == tokens[2][0]:
          maxsplit = 6 #decoys name.txt file
        elif 12 == len(tokens[0]) and "T"==tokens[0][0] and 9 != len(tokens[1]):
          maxsplit = 2 #dbgen file
        elif tokens[0] == 'name.txt': #new dbstart
          maxsplit = 6
        elif tokens[0] == 'name.cxcalc.txt': #new dbstart cxcalc version
          maxsplit = 8
        else:
          maxsplit = 4 #ligands name.txt file
        #do the split again, with maxsplit
        tokens = firstLine.split(None, maxsplit) #None forces use of whitespace
        if tokens[0] == 'name.txt': #new dbstart
          self.name = tokens[2]
          self.protName = "none"
          #self.smiles = string.strip(tokens[3])
          #self.longname = string.strip(tokens[5])
          self.smiles = tokens[3].strip()
          self.longname = tokens[5].strip()
        elif tokens[0] == 'name.cxcalc.txt': #new dbstart
          self.name = tokens[2]
          self.protName = tokens[3]
          self.smiles = tokens[4].strip()
          self.longname = tokens[7].strip()
        elif 7 == len(tokens): #decoys name.txt file
          self.name = tokens[1]
          self.protName = tokens[2]
          self.smiles = tokens[3].strip()
          self.longname = tokens[6].strip()
        elif 5 == len(tokens): #ligands name.txt file
          self.name = tokens[1]
          self.protName = "none"
          self.smiles = tokens[2].strip()
          self.longname = tokens[4].strip()
        elif 3 == len(tokens): #dbgen name.txt file
          self.name = tokens[0]
          self.protName = "none"
          self.smiles = tokens[1].strip()
          self.longname = tokens[2].strip()
        #print self.name, self.protName, self.smiles, self.longname #debug
      except StopIteration: #end of file
        namefile.close()
    if mol2fileName is not None:
      if mol2fileName.endswith(".gz"):
        mol2file = gzip.GzipFile(mol2fileName, 'r')
      else:
        mol2file = open(mol2fileName, 'r')
      self.phase = 0
      try:
        for line in mol2file:
          self.processLine(line)
      except StopIteration:
        mol2file.close()
      finally:
        self.xyzCount += 1 #really this needs to be done
        self.origXyzCount = self.xyzCount
        while len(self.inputEnergy) < self.xyzCount:
          self.inputEnergy.append(9999.99)
        while len(self.inputHydrogens) < self.xyzCount:
          self.inputHydrogens.append(0)
    #print self.atomName,    self.atomType, self.atomCharge
    #print self.bondStart,     self.bondEnd, self.bondType
    #print self.xyzCount

  def processLine(self, line):
    '''reads a single line, processes it'''
    if line[:17] == "@<TRIPOS>MOLECULE":
      self.phase = 1 #header phase
    elif line[:13] == "@<TRIPOS>ATOM":
      self.phase = 2 #atoms phase
      self.xyzCount += 1 #first one is numbered 0
      self.atomXyz.append([]) #list of lists
    elif line[:13] == "@<TRIPOS>BOND":
      self.phase = 3 #bonds phase
    elif line[:9] == "@<TRIPOS>":
      self.phase = 0 #fake phase that reads nothing...
    elif line[0] == '#': #comment line, ignore
      pass
    elif len(line) == 1: #comment line, ignore
      pass
    else:
      if 1 == self.phase:
        if self.name == "fake":
          self.name = line.strip()
        #tokens = string.split(string.strip(line))
        tokens = line.strip().split()
        if len(line) > 1 and tokens[0][0:7] == "mmff94s":
          self.inputEnergy.append(float(tokens[2]))
          self.inputHydrogens.append(0)
          self.phase = 0 #rest of header ignored
      elif 2 == self.phase:
        #tokens = string.split(line)
        tokens = line.strip().split()
        if 0 == self.xyzCount: #only read in some stuff for first mol2
          self.atomNum.append(int(tokens[0]))
          self.atomName.append(tokens[1])
          self.atomType.append(tokens[5])
          self.atomCharge.append(float(tokens[-1])) #last column is charge
          self.atomBonds.append([]) #start new list
        self.atomXyz[self.xyzCount].append( \
           (float(tokens[2]), float(tokens[3]), float(tokens[4])))
      elif 3 == self.phase and 0 == self.xyzCount:
        #bonds only read in for first molecule, assumed the same after
        #tokens = string.split(line)
        tokens = line.strip().split()
        self.bondNum.append(int(tokens[0]))
        self.bondStart.append(int(tokens[1]))
        self.bondEnd.append(int(tokens[2]))
        self.bondType.append(tokens[3])
        self.atomBonds[int(tokens[1]) - 1].append( \
                       (int(tokens[2]) - 1, tokens[3]))
        self.atomBonds[int(tokens[2]) - 1].append( \
                       (int(tokens[1]) - 1, tokens[3]))

  def copy(self):
    '''returns new mol2 object'''
    newM = Mol2()
    newM.name = self.name
    newM.protName = self.protName
    newM.atomNum = self.atomNum
    newM.atomName = self.atomName
    newM.atomXyz = self.atomXyz[:] #only this
    newM.inputEnergy = self.inputEnergy[:] #and this
    newM.xyzCount = self.xyzCount #and this
    newM.origXyzCount = self.origXyzCount #and this
    newM.inputHydrogens = self.inputHydrogens[:] #and this are to be messed with
    newM.atomType = self.atomType
    newM.atomCharge = self.atomCharge
    newM.atomBonds = self.atomBonds
    newM.bondNum = self.bondNum
    newM.bondStart = self.bondStart
    newM.bondEnd = self.bondEnd
    newM.bondType = self.bondType
    newM.smiles = self.smiles
    newM.longname = self.longname
    try:
      newM.colorConverter = self.colorConverter
      newM.dockNum = self.dockNum
      newM.colorNum = self.colorNum
    except AttributeError:
      pass #this is fine
    newM.bondDists = None #don't copy this, regenerate
    return newM #new Mol2 object that can be manipulated

  def initFromDb2Lines(self, mlines):
    '''reads data from lines in the db2 file. start with a blank, init obj.
    mlines is lines starting with M'''
    tokens = []
    for mline in mlines[0:4]: #1st 4 mlines only
      tokens.append(string.split(mline))
    self.name = tokens[0][1] #2d array, 0th line, 1st token
    self.protName = tokens[0][2]
    self.smiles = tokens[2][1]
    self.longname = tokens[3][1]
    self.atomNum = []
    self.atomName = [] #all this stuff might need to get optimized if lists suck
    self.atomXyz = [] #kept for every conformation
    self.inputEnergy = [] #kept for every conformation
    self.inputHydrogens = [] #for every conf. 0 means input. 1 means reset.
                            #2 means rotated. 3 means mix-n-match so who knows
    self.atomType = []
    self.atomCharge = [] #read in but overriden by solv data
    self.atomBonds = []
    self.bondNum = []
    self.bondStart = []
    self.bondEnd = []
    self.bondType = []
    self.xyzCount = -1
    self.origXyzCount = -1

  def keepConfsOnly(self, first, last):
    ''' delete input confs outside of the first, last range'''
    if last > self.xyzCount: #can't copy beyond end
      last = self.xyzCount
    self.xyzCount = last - first
    self.atomXyz = self.atomXyz[first:last]
    self.inputEnergy = self.inputEnergy[first:last]
    self.inputHydrogens = self.inputHydrogens[first:last]

  def bondsBetween(self, atomNum, atomOther):
    '''returns number of bonds between any 2 atoms (see bondsBetweenActual)
    atomNum and otherNum is the mol2 numbering (1-based). '''
    actualNum = atomNum - 1 #assume 1-based to 0-based conversion
    actualOther = atomOther - 1  #assume 1-based to 0-based conversion
    return self.bondedToActual(actualNum, atomOther)

  def bondsBetweenActual(self, actualNum, actualOther):
    '''returns the number of bonds between the two atom numbers specified.
    directly bonded => 1
    one atom in between => 2
    two atoms => 3, etc.'''
    if self.bondDists is None: #generate them
      self.calcBondDists()
    row = self.bondDistsOrderKeys[actualNum]
    col = self.bondDistsOrderKeys[actualOther]
    return self.bondDists[row][col]

  def distFromAtoms(self, atoms):
    '''using a particular set of atoms (usually rigid component) as 0,
    find the bond distance to all other atoms'''
    neighbors = {}
    for atomNum in range(len(self.atomNum)):
      neighbors[atomNum] = []
      for otherBond in self.atomBonds[atomNum]:
        neighbors[atomNum].append((otherBond[0], 1))
    dists = shortestpaths.shortestPaths(neighbors.keys(), neighbors, 0, atoms)
    return dists

  def calcBondDists(self):
    '''uses floyd warshall all pairs shortest paths algorithm to generate
    the # of bonds between all pairs of atoms. cache results and don't redo'''
    neighbors = {}
    for atomNum in range(len(self.atomNum)):
      neighbors[atomNum] = []
      for otherBond in self.atomBonds[atomNum]:
        neighbors[atomNum].append((otherBond[0], 1))
    distances, orderKeys = floydwarshall.floydWarshall(neighbors)
    self.bondDists = distances
    self.bondDistsOrderKeys = orderKeys
    #no return, we're done here

  def bondedTo(self, atomNum, firstName, bondsAway=1, lastBond=None, \
               returnAtom=False):
    '''returns true if atomNum has a bond to any other atom with a name starting
    with firstName. atomNum is the mol2 numbering (1-based). bondsAway
    controls how far the search proceeds before checking. must be exact.'''
    actualNum = atomNum - 1 #assume 1-based to 0-based conversion
    return self.bondedToActual(actualNum, firstName, bondsAway, lastBond, \
                               returnAtom)

  def bondedToAll(self, atomNum, firstName, bondsAway=1, lastBond=None):
    '''returns all atoms bonded, for other functions to process. calls
    bondedToActualAll, see for more documentation'''
    actualNum = atomNum - 1 #assume 1-based to 0-based conversion
    return self.bondedToActualAll(actualNum, firstName, bondsAway, lastBond)

  def bondedToActual(self, actualNum, firstName, bondsAway=1, lastBond=None, \
                     returnAtom=False):
    '''returns true if actualNum has a bond to any other atom with a name
    starting with firstName. actualNum is the 0-based numbering. bondsAway
    controls how far the search proceeds before checking. must be exact.
    lastBond controls the type of the lastBond found (to the one just before
     the ones being checked) useful for finding only certain types of bonds.
    if returnAtom is True, it returns the atom number of the atom that matched.
    '''
    bondedAwayNums = self.bondedToActualAll(actualNum, firstName, bondsAway, \
                                            lastBond)
    for anAtomNum in bondedAwayNums[bondsAway]:
      #if -1 != string.find(self.atomType[anAtomNum], firstName):
      if -1 != self.atomType[anAtomNum].find(firstName):
        if not returnAtom:
          return True #found a matching atom at the other end of a bond
        else: #returnAtom is True
          return True, self.atomNum[anAtomNum]
    if not returnAtom:
      return False #no bonded atoms with that first letter
    else: #returnAtom is True
      return False, False

  def bondedToActualAll(self, actualNum, firstName, bondsAway=1, \
                        lastBond=None):
    '''returns all atoms a certain number of bonds away, obeying lastBond as
    bondedToActual documents'''
    bondedAwayNums = collections.defaultdict(list) #defaults to empty list
    bondedAwayNums[0].append(actualNum)
    checked = 0
    while checked < bondsAway: #check up to bondsaway bonds away from start
      for startNum in bondedAwayNums[checked]: #for each atom in cur level
        for otherBond in self.atomBonds[startNum]:
          #check to make sure otherBond[0] not in any previous list
          okayToAdd = True
          #for checkList in bondedAwayNums.itervalues():
          for checkList in iter(bondedAwayNums.values()):
            if otherBond[0] in checkList:
              okayToAdd = False
          if okayToAdd and (lastBond is not None) and (checked +1 == bondsAway):
            if lastBond != "*": #same as none, basically allow anything
              if -1 == otherBond[1].find(lastBond): #-1 means no match
                okayToAdd = False #which means this bond isn't correct
          if okayToAdd:
            bondedAwayNums[checked + 1].append(otherBond[0])
      checked += 1 #move 1 more away
    return bondedAwayNums

  def isAtomBondedOtherThan(self, atomNum, count, otherThan):
    '''for each atom, if it is bonded to any of [count] atoms that are not of
    type [otherThan], return true'''
    bondedAwayNums = self.bondedToAll(atomNum, "")[1] #only want 1 bond away
    otherThanCount = 0
    for atomNumActual in bondedAwayNums:
      if self.atomType[atomNumActual] not in otherThan:
        otherThanCount += 1
    if otherThanCount not in count:
      return True
    else:
      return False

  def convertDockTypes(self, parameterFileName=None):
    '''adds self.dockNum to each atom record based on sybyl2dock'''
    dockConverter = sybyl2dock.AtomConverter(parameterFileName)
    self.dockNum = [] #new data on each dock atom number
    for atomNumber in self.atomNum:
      self.dockNum.append(dockConverter.convertMol2atomNum(self, atomNumber))

  def addColors(self, parameterFileName=None):
    '''adds self.colorNum to each atom record based on rules'''
    colorConverter = atom_color_table.ColorTable(parameterFileName)
    self.colorConverter = colorConverter #save for later use in output
    self.colorNum = [] #map from atom to colors
    for atomNumber in self.atomNum:
      self.colorNum.append(colorConverter.convertMol2color(self, atomNumber))

  def countConfs(self):
    '''returns the number of conformations in the file'''
    return len(self.atomXyz)

  def getXyz(self, xyzCount, atomNum):
    '''returns the xyz for an atom number and an xyz count (conformation)'''
    atomIndex = self.atomNum.index(atomNum)
    return self.atomXyz[xyzCount][atomIndex]

  def getXyzManyConfs(self, xyzCounts, atomIndex):
    '''returns a lits of xyzs for many confs of one atom'''
    xyzConfs = []
    for xyzCount in xyzCounts:
      xyzConfs.append(self.atomXyz[xyzCount][atomIndex])
    return xyzConfs

  def getRMSD(self, xyzOne, xyzTwo):
    '''calculates just the rmsd of the two conformations'''
    sumSquared = 0.0
    for atomIndex in range(len(self.atomXyz[xyzOne])):
      sumSquared += geometry.distL2Squared3( \
                        self.atomXyz[xyzOne][atomIndex],
                        self.atomXyz[xyzTwo][atomIndex])
    rmsd = (sumSquared / len(self.atomXyz[xyzOne])) ** 0.5
    return rmsd

  def getRMSDtable(self, forceRedo=False):
    '''for each conformation (xyzcount) to all others, find the rmsd. return
    as dictionary of dictionaries of rmsds'''
    if self.rmsdTable is None or forceRedo:
      self.rmsdTable = {}
      self.rmsdList = []
      for xyzCount in range(len(self.atomXyz)):
        self.rmsdTable[xyzCount] = {} #set to empty dicts first
      for xyzCount in range(len(self.atomXyz)):
        for otherXyz in range(xyzCount + 1, len(self.atomXyz)): #just half
          rmsd = self.getRMSD(xyzCount, otherXyz)
          self.rmsdTable[xyzCount][otherXyz] = rmsd
          self.rmsdTable[otherXyz][xyzCount] = rmsd
          self.rmsdList.append((rmsd, otherXyz, xyzCount))
      self.rmsdList.sort(key=operator.itemgetter(0))
    return self.rmsdTable

  def getRMSDlist(self):
    '''gets the rmsds between all pairs as a list of rmsd, conf, conf tuples'''
    self.getRMSDtable()
    return self.rmsdList

  def getRMSDclusters(self, rmsdCutoff=None, numClusters=1):
    '''uses the rmsdlist to make clusters of conformations based on rmsd.
    goes until either the rmsdCutoff is reached or numClusters is reached.
    using the numClusters will make this run very slowly.
    uses single linkage to make a new cluster.'''
    self.getRMSDtable() #make the table, or ensure it is made
    #self.rmsdList is a tuple of (rmsd, conf, conf)
    clusters = unionfind2.unionFind()
    for xyzCount in range(len(self.atomXyz)):
      clusters.find(xyzCount) # initialize all these to singleton clusters
    if rmsdCutoff is None:
      rmsdCutoff = self.rmsdList[-1][0] + 1.0 #make it never happen
    for rmsdTuple in self.rmsdList:
      if rmsdTuple[0] > rmsdCutoff:
        break #quit combining things!
      clusters.union(rmsdTuple[1], rmsdTuple[2])
    return clusters.toLists()

  def getRMSDclustersAll(self, rmsdCutoff=None, numClusters=1):
    '''uses the rmsdlist to make clusters of conformations based on rmsd.
    goes until either the rmsdCutoff is reached or numClusters is reached.
    using the numClusters will make this run very slowly. uses ALL linkage not
    single linkage.'''
    self.getRMSDtable() #make the table, or ensure it is made
    #self.rmsdList is a tuple of (rmsd, conf, conf)
    #self.rmsdTable is a dict of [conf][conf] -> rmsd
    clusters = unionfind2.unionFind()
    for xyzCount in range(len(self.atomXyz)):
      clusters.find(xyzCount) #init
    if rmsdCutoff is None:
      rmsdCutoff = self.rmsdList[-1][0] + 1.0 #make it never happen
    for rmsdTuple in self.rmsdList:
      if rmsdTuple[0] > rmsdCutoff:
        break #quit combining things!
      #have to do all linkage not just single.. oh my
      if clusters.different(rmsdTuple[1], rmsdTuple[2]): #otherwise already join
        combine = True
        clusterOne = clusters.getList(rmsdTuple[1])
        clusterTwo = clusters.getList(rmsdTuple[2])
        #print clusterOne, clusterTwo,
        for clusterOneRep in clusterOne:
          for clusterTwoRep in clusterTwo:
            thisRMSD = self.rmsdTable[clusterOneRep][clusterTwoRep]
            #print thisRMSD,
            if thisRMSD > rmsdTuple[0]: #means we can't combine yet
              combine = False
              break
          if not combine:
            break
        #print combine
        if combine:
          clusters.union(rmsdTuple[1], rmsdTuple[2])
    return clusters.toLists()

  def divisiveClustering(self):
    '''takes all conformations. bisects them along the longest dimension
    (in N*atoms*3 space). Repeat on the biggest remaining cluster until there
    are numClusters left, return the clusters as lists.'''
    numClusters = min(20, int(self.origXyzCount/3.))
    #print numClusters #debugging, find out target # of clusters
    clusters = divisive_clustering.divisiveClustering(self.atomXyz, numClusters)
    return clusters

  def addSolvDataPartialCharges(self, partialCharges):
    '''adds partial charge data from list, must be in correct order'''
    for count in range(len(self.atomCharge)):
      self.atomCharge[count] = partialCharges[count]
    #that's it. no return.

  def writeMol2File(self, outFile, whichXyz=None):
    '''writes the data to an already open file. don't close it.'''
    if whichXyz is None:
      whichXyz = range(self.xyzCount)
    for oneXyz in whichXyz:
      outFile.write("@<TRIPOS>MOLECULE\n")
      if self.protName == "fake": #don't write fake
        outFile.write(self.name + "\n")
      else:
        outFile.write(self.name + " " + self.protName + "\n")
      outFile.write("%5d %5d %5d %5d %5d\n" % (len(self.atomNum), \
                    len(self.bondNum), 0, 0, 0))
      outFile.write("SMALL\nUSER_CHARGES\n\n")
      outFile.write("mmff94s_NoEstat = %5.2f\n" % self.inputEnergy[oneXyz])
      outFile.write("@<TRIPOS>ATOM\n")
      for oneAtom in range(len(self.atomNum)):
        outFile.write( \
          "%7d %6s    % 8.4f  % 8.4f  % 8.4f %5s     1 <0>       % 8.4f\n" % \
          (self.atomNum[oneAtom], string.ljust(self.atomName[oneAtom], 6), \
           self.atomXyz[oneXyz][oneAtom][0], \
           self.atomXyz[oneXyz][oneAtom][1], \
           self.atomXyz[oneXyz][oneAtom][2], \
           string.ljust(self.atomType[oneAtom], 5), \
           self.atomCharge[oneAtom]))
      outFile.write("@<TRIPOS>BOND\n")
      for oneBond in range(len(self.bondNum)):
        outFile.write( \
          "%6d%5d%5d %2s\n" % \
          (self.bondNum[oneBond], self.bondStart[oneBond], \
           self.bondEnd[oneBond], string.ljust(self.bondType[oneBond], 2)))

  def writeMol2(self, outName, whichXyz=None):
    '''writes the data to a mol2 file.'''
    outFile = open(outName, 'w')
    self.writeMol2File(outFile, whichXyz)
    outFile.close()

  def deleteBond (self, bondInd):
    '''deletes the bond at index bondInd and adjust the rest of the indexing'''
    nbonds = len(self.bondNum)
    bondi = bondInd+1;
    while (bondi < nbonds):
      self.bondNum[bondi] = self.bondNum[bondi]-1
      bondi = bondi+1

    del self.bondNum[bondInd]
    del self.bondStart[bondInd]
    del self.bondEnd[bondInd]
    del self.bondType[bondInd]

  def deleteAtom (self, atomInd):
    '''deletes atom number atomInd from all tables and adjust the rest of the 
    indexing and bond numbers accordingly'''

    #removes atomnum and adjusts numbering
    #print atomInd
    natoms = len(self.atomNum)
    atomn = self.atomNum[atomInd]
    atomi = atomInd+1;
    while (atomi < natoms):
      self.atomNum[atomi] = self.atomNum[atomi]-1
      atomi = atomi+1

    #removes atom from atoms tables
    del self.atomNum[atomInd]
    del self.atomType[atomInd]
    del self.atomBonds[atomInd]
    del self.atomName[atomInd]
    del self.atomCharge[atomInd]
    del self.dockNum[atomInd]
    del self.colorNum[atomInd]
    for x in self.atomXyz:
      del x[atomInd]
 
    #delete any bond containing the involved atom
    nbonds = len(self.bondNum)
    bondi = 0
    while (bondi < nbonds):
      if (self.bondStart[bondi]==atomn or self.bondEnd[bondi]==atomn):
        self.deleteBond(bondi)
        bondi = bondi-1
        nbonds = len(self.bondNum)
      bondi = bondi+1

    #adjust bonds atom numbering
    nbonds = len(self.bondNum)
    bondi = 0
    while (bondi < nbonds):
      if (self.bondStart[bondi] >= atomn):
        self.bondStart[bondi] = self.bondStart[bondi] - 1
      if (self.bondEnd[bondi] >= atomn):
        self.bondEnd[bondi] = self.bondEnd[bondi] - 1
      bondi = bondi+1

    natoms = len(self.atomNum)
    atomi = 0
    while (atomi < natoms):
      nbonds = len(self.atomBonds[atomi])
      bondi = 0
      while (bondi < nbonds):
        if (self.atomBonds[atomi][bondi][0] >= atomInd):
          tmpBond = (self.atomBonds[atomi][bondi][0]-1, self.atomBonds[atomi][bondi][1])
          self.atomBonds[atomi][bondi]=tmpBond
        bondi = bondi+1
      atomi = atomi+1

  def recolorCovalentAttachment (self, covAtomType):
    '''assumes the covalent dummy atom is already removed, and the covalent attachmet
    atom type is already set to 'cov' will color it "8" and it's bonded neighbours: 
    9, 10, 11'''
    
    #find the covalent attachment atom
    for atom in self.atomType:
      if (atom == "cov"):
        cov_index = self.atomType.index(atom)

    color = 8 #well known covalent color :)
    self.colorNum[cov_index] = color 
    for bond in self.atomBonds[cov_index]:
      color = color + 1
      neighbor = bond[0]
      self.colorNum[neighbor] = color
    #reset the original covalent atom type returned by removeCovalentDummyAtom
    self.atomType[cov_index] = covAtomType

  def removeCovalentDummyAtom (self):
    '''if this is a covalent ligand marked by a SiH3 dummy atom, removes the Si
    emoves the 3H it is connected to, their bonds, and colors the covalent
    attachment atom approprietly'''

    indicesList = []

    #find out what the Si atom index is and mark it for deletion
    for atom in self.atomType:
      if (atom == "Si"):
        si_index = self.atomType.index(atom)
        self.atomType[si_index] = "del"

    #find out the indices of it's hydrogen neighbours and cov attach point
    
    for bond in self.atomBonds[si_index]:
      if self.atomType[bond[0]] == "H":
        self.atomType[bond[0]] = "del"
        print (str(bond[0]) + " to be deleted")
      else:
        covAtomType = self.atomType[bond[0]] 
        self.atomType[bond[0]] = "cov"
        print (str(bond[0]) + " is the cov attach")
        # explicitly remove the bond to the covalent attachment point    
        nbonds = len(self.atomBonds[bond[0]])
        bondi = 0
        while (bondi < nbonds):
          print (nbonds)
          print (bond[0])
          print (bondi)
          tmpbond = self.atomBonds[bond[0]][bondi]
          if (tmpbond[0] == si_index):
            del self.atomBonds[bond[0]][bondi]
            nbonds = nbonds - 1
            bondi = bondi - 1 
          bondi  = bondi + 1

    #find the indices of the atoms to be removed, so we can 
    #also remove them from the output.solv output
    natoms = len(self.atomNum)  
          
    #now actually go ahead and delete these atoms
    atomindex = 0
    while (atomindex < natoms):
      if self.atomType[atomindex] == "del":
        indicesList.append(atomindex)
        self.deleteAtom(atomindex)
        atomindex = atomindex-1
        natoms = len(self.atomNum)
        print ("natoms " + str(natoms))
      atomindex = atomindex+1      

    return covAtomType, indicesList
  
def readDockMol2file(mol2fileName, recdes=False, ligdes=False, charge=False, \
    elec=False):
  '''reads a dock output mol2 file, since each ligand has different connectivity
  this returns a list of Mol2 data classes instead of just one big one.
  if desired, can return the receptor desolvation as a list and/or
  the polar ligand desolvation scores as well'''
  mol2lines = []
  mol2data = []
  mol2rd = []
  mol2ld = []
  mol2charge = []
  mol2elec = []
  mol2file = open(mol2fileName,	'r')
  mol2names = []
  for mol2line in mol2file:
    if mol2line[:17] == "@<TRIPOS>MOLECULE":
      if len(mol2lines)	> 0:
       	mol2data.append(Mol2(mol2text=mol2lines))
      mol2lines	= []
    if mol2line[0] != "#" and len(mol2line) > 1:
      mol2lines.append(mol2line)
    if mol2line[:32] == '##########                 Name:':
      mol2names.append(string.split(mol2line)[2])
    if mol2line[:32] == '##########  Ligand Polar Desolv:':
      mol2ld.append(float(string.split(mol2line)[4]))
    if mol2line[:32] == '########## Receptor Desolvation:':
      mol2rd.append(float(string.split(mol2line)[3]))
    if mol2line[:32] == '##########        Ligand Charge:':
      mol2charge.append(float(string.split(mol2line)[3]))
    if mol2line[:32] == '##########        Electrostatic:':
      mol2elec.append(float(string.split(mol2line)[2]))
  if len(mol2lines) > 0:
    mol2data.append(Mol2(mol2text=mol2lines))
  for count, oneMol2 in enumerate(mol2data):
    oneMol2.name = mol2names[count]
  retList = [mol2data]
  if recdes:
    retList.append(mol2rd)
  if ligdes:
    retList.append(mol2ld)
  if charge:
    retList.append(mol2charge)
  if elec:
    retList.append(mol2elec)
  return tuple(retList)

#if 0 == string.find(sys.argv[0], "mol2.py"):
if 0 == sys.argv[0].find("mol2.py"):
  #if called from commandline, read in all arguments as mol2 files and print
  #out statistics about the molecule or confirm that it was read in who knows
  for mol2file in sys.argv[1:]:
    mol2data = Mol2(mol2file)
    print (mol2data.name)
    print (len(mol2data.atomName), len(mol2data.bondStart)) #debugging???
    mol2data.convertDockTypes()
    print (mol2data.dockNum)
