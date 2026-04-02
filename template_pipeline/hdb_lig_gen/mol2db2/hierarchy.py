#!/usr/bin/env python2.7

#Ryan G. Coleman
#uses mol2 file to generate a hierarchy

import string
import sys
from unionfind2 import unionFind
import geometry
import buckets
import gzip
import operator
import math
import time
import shortestpaths

def printClusterHelper(clusterList):
  '''stupid function used for debugging, prints list of pymol out.???.mol2 lines
  to copy/paste and run to see what the clusters are.
  run mol2hydroxyls.py -r and mol2tomultimol2.py first to get out.???.mol2 files
  '''
  for clusters in clusterList:
    strout = "pymol "
    for conf in clusters:
      strout = stout + "out." + string.zfill(conf, 3) + ".mol2 "
    print (strout)

def computeBreaks(limitError, options):
  '''3 diff requirements, make sure we break it into enough pieces to meet them
  all.'''
  #have to break the atomXyz into multiple sets so the hierarchy isn't too big
  try:
    breaksS = int(math.ceil(limitError.getSets() / float(options.limitset)))
  except TypeError:  # means None was used
    breaksS = 1
  try:
    breaksC = int(math.ceil(limitError.getConfs() / float(options.limitconf)))
  except TypeError:  # means None was used
    breaksC = 1
  try:
    breaksX = int(math.ceil(limitError.getCoords() / float(options.limitcoord)))
  except TypeError:  # means None was used
    breaksX = 1
  #print breaksS, breaksC, breaksX  # see which breaks is higher
  breaks = max(breaksS, breaksC, breaksX)  # use the max of any of these
  return breaks

class TooBigError(Exception):
  '''error raised when the hierarchy has too many conformations of input
  after the hydroxyls have been rotated.'''

  def __init__(self, confs, sets, coords):
    self.confs = confs
    self.coords = coords
    self.sets = sets

  def __str__(self):
    return repr(self.confs) + ", " + repr(self.coords) + ", " + repr(self.sets)

  def getConfs(self):
    '''actually used to figure out how many sub-groups to split input confs'''
    return self.confs

  def getCoords(self):
    '''actually used to figure out how many sub-groups to split input confs'''
    return self.coords

  def getSets(self):
    '''actually used to figure out how many sub-groups to split input confs'''
    return self.sets

class Hierarchy(object):
  '''uses data from a mol2 file to make a hierarchy of conformations.
  the following constants are used when writing out the confs/groups and are
  based on the 80 character limit in fortran. yeah seriously.
  they might change if something serious happens but it is better that they
  are here than hardcoded several times later
  these are floats so that the division works'''
  grGrPerLine = 17.  # group -> group children per line in output
  grCoPerLine = 9.  # group -> conf
  coCoPerLine = 9.  # conf -> conf
  coSePerLine = 8.  # conf -> set

  def __init__(
      self, mol2data, clashDecider, tolerance=0.001, verbose=False,
      timeit=False, limitset=9999999999, limitconf=9999999999,
      limitcoord=9999999999, solvdata=None):
    '''takes a mol2data class as input. makes a hierarchy.'''
    if solvdata is not None:
      self.solvdata = solvdata
    if timeit:
      startTime = time.time()
    #first step is to count the number of positions each atom has.
    #the tolerance is taken into account here and only here.
    totalCoords = len(mol2data.atomXyz) * len(mol2data.atomXyz[0])
    if verbose:
      print ("total number of sets (complete confs):", len(mol2data.atomXyz))
    if len(mol2data.atomXyz) > limitset:  # quit now, way too many sets
      if verbose:
          print ("Total number of sets too high %d > %d" % (len(mol2data.atomXyz), limitset))
      raise TooBigError(None, len(mol2data.atomXyz), totalCoords)
    if len(mol2data.atomXyz) > 50:
      if verbose:
        print ("using faster count positions algorithm for large data")
      self._countPositions(mol2data.atomXyz, tolerance, verbose)
    else:
      if verbose:
        print ("using default count positions algorithm for smaller data")
      self._countPositionsFewPoints(mol2data.atomXyz, tolerance)
    if timeit:
      countTime = time.time()
      print ("time to count unique positions:", countTime-startTime)
    if verbose:
      print ("unique positions, atoms:", self.posCount, len(mol2data.atomXyz))
    if totalCoords > limitcoord:
      if verbose:
          print ("Total number of coords too high %d > %d" % (totalCoords, limitcoord))
      raise TooBigError(None, len(mol2data.atomXyz), totalCoords)
      #this breaks out of the init stage, needs fewer confs to be passed in.
    #the rigid component is the biggest set of bonded non-moving atoms
    self._findRigidComponent(mol2data.atomBonds)  # also uses self.posCount
    if timeit:
      rigidTime = time.time()
      print ("time to find rigid component:", rigidTime-countTime)
    if verbose:
      print ("rigid atoms, others:", self.rigidComponent, self.atomsNotAssigned)
    #new algorithm, find bonded atoms that move together, put in conformations
    self._findRigidHeavy(mol2data.atomType)
    self.heavyAtomNums = None
    self._setHeavy(mol2data.atomType)
    self._findConformations(mol2data.atomBonds, mol2data.atomXyz)
    self._findSets()  # puts conformations in sets
    if timeit:
      flexTime = time.time()
      print ("time to find flexible components:", flexTime-rigidTime)
    if verbose:
      print ("total number of confs:", self.confNums[-1])
    if self.confNums[-1] > limitconf:
      raise TooBigError(
          self.confNums[-1], len(mol2data.atomXyz), totalCoords)
      #this breaks out of the init stage, needs fewer confs to be passed in.
    #now want to actually put atom positions into hierarchy groups
    self._assignCoords(mol2data.atomXyz)
    if timeit:
      assignCoordsTime = time.time()
      print ("time to assign coords:", assignCoordsTime-flexTime)
    self._identifyClashSetnums(clashDecider, mol2data)
    if timeit:
      afterClash = time.time()
      print ("time to identify clash sets:", afterClash-assignCoordsTime)
    if verbose:
      print ("number of broken/clashed sets:", len(self.brokenSets))
    #the mol2data is needed during output so save it.
    self.mol2data = mol2data
    if timeit:
      afterXyz = time.time()
      print ("time to identify conf atoms:", afterXyz - afterClash)
    self.clusters = None  # used to detect if clustering/clouding was done
    self._makeClouds()  # highest level of ligand sampling
    if timeit:
      afterClouds = time.time()
      print ("time to make clouds:", afterClouds - afterXyz)

  def _countPositions(self, xyzData, tolerance, verbose=False):
    '''for a list of list of xyz data, count the number of positions each
    atom takes based on the tolerance and the distance. tolerance is compared
    to the euclidean difference squared to determine if a position is equal.
    actually uses a clustering algorithm and uses a unionfind data structure.'''
    self.posCount = []
    self.posClusters = []  # just save all the data since we made it
    self.posClusterLists = []  # just save all the data since we made it
    tolerance2 = tolerance ** 2.  # square the tolerance since it is compared
    for oneSet in range(len(xyzData[0])):  # goes from 0 to atom count
      #if verbose:
      #  print oneSet, " atom positions being calculated"
      clusters = unionFind()
      xyzList = []
      for oneIndex in range(len(xyzData)):  # 0 to number of positions (mol2#s)
        clusters.find(oneIndex)  # initiate each position
        xyzList.append(xyzData[oneIndex][oneSet])
      bucket = buckets.Bucket3d(xyzList, tolerance)  # constructor to make fast
      bucket.getWithinCluster(clusters)
      #for pointA, pointB in bucket.getWithin(clusters):
      #  clusters.union(pointA, pointB)
      tempLists = clusters.toLists()
      self.posCount.append(len(tempLists))
      self.posClusters.append(clusters)
      self.posClusterLists.append(tempLists)

  def _countPositionsFewPoints(self, xyzData, tolerance):
    '''for a list of list of xyz data, count the number of positions each
    atom takes based on the tolerance and the distance. tolerance is compared
    to the euclidean difference squared to determine if a position is equal.
    actually uses a clustering algorithm and uses a unionfind data structure.'''
    self.posCount = []
    self.posClusters = []  # just save all the data since we made it
    self.posClusterLists = []  # just save all the data since we made it
    tolerance2 = tolerance ** 2.  # square the tolerance since it is compared
    for oneSet in range(len(xyzData[0])):  # goes from 0 to atom count
      clusters = unionFind()
      xyzList = []
      for oneIndex in range(len(xyzData)):  # 0 to number of positions (mol2#s)
        clusters.find(oneIndex)  # initiate each position
        xyzList.append(xyzData[oneIndex][oneSet])
      for oneIndex in range(len(xyzData)):  # 0 to positions
        oneXyz = xyzList[oneIndex]
        for twoIndex in range(oneIndex+1, len(xyzData)):
          # count from oneIndex to positions
          if geometry.distL2Squared3(oneXyz, xyzList[twoIndex]) < tolerance2:
            clusters.union(oneIndex, twoIndex)
      tempLists = clusters.toLists()
      self.posCount.append(len(tempLists))
      self.posClusters.append(clusters)
      self.posClusterLists.append(tempLists)

  def _findRigidComponent(self, atomBonds):
    '''uses bond and position count information to find largest set of atoms
    that don't move. this is the rigid component. set into self.rigidComponent
    also find the complement of atomnums and the rigid component and set into
    self.atomsNotAssigned for use later'''
    clusters = unionFind()
    for atomNum in range(len(self.posCount)):
      if 1 == self.posCount[atomNum]:
        for otherNum, bondType in atomBonds[atomNum]:
          if 1 == self.posCount[otherNum]:
            clusters.union(atomNum, otherNum)
    maxSize = 0
    maxCluster = None
    clusterLists = clusters.toLists()
    for clusterList in clusterLists:
      if len(clusterList) > maxSize:
        maxSize = len(clusterList)
        maxCluster = clusterList
    self.rigidComponent = maxCluster
    self.atomsAssigned = set(self.rigidComponent)
    self.atomsNotAssigned = set()
    for atomNum in range(len(self.posCount)):
      if atomNum not in self.rigidComponent:
        self.atomsNotAssigned.add(atomNum)

  def _findRigidHeavy(self, atomTypes):
    '''counts the heavy atoms in the rigid component and puts in
    self.heavyRigidCount'''
    self.heavyRigidCount = 0
    self.heavyRigidAtomNums = []
    for atomNum in self.atomsAssigned:
      if atomTypes[atomNum].find('H') == -1:
        self.heavyRigidAtomNums.append(atomNum)
        self.heavyRigidCount += 1
    #print self.heavyRigidCount

  def _setHeavy(self, atomTypes):
    '''for all atoms, finds the heavy ones, put in self.heavyAtomNums, return'''
    if self.heavyAtomNums is None:  # only do this once, it never changes
      self.heavyAtomNums = []
      for atomNum in range(len(atomTypes)):
        if atomTypes[atomNum].find('H') == -1:
          self.heavyAtomNums.append(atomNum)
    return self.heavyAtomNums

  def _findConformations(self, atomBonds, xyzData):
    '''uses bond and xyzs to figure out what sets of neighboring atoms move
    together and assign them to conformations and assign each set a specific
    bunch of conformations.
    self.rigidComponent is the list of atom numbers for the rigid comp
    self.atomsAssigned is the set of atom numbers for the rigid comp (@start)
    self.atomsNotAssigned is the rest of the atom numbers'''
    self.confNums = [1]  # rigid starts
    self.confAtoms = {}  # maps to atom numbers
    self.confAtoms[1] = list(self.atomsAssigned)
    self.confInput = {}  # maps to the input xyz lists
    self.confInput[1] = range(len(xyzData))
    confClusters = {}
    for atomNum in self.atomsNotAssigned:
      for listInputs in self.posClusterLists[atomNum]:
        tupleInputs = tuple(listInputs)  # can't use lists as keys
        if tupleInputs not in confClusters.keys():
          confClusters[tupleInputs] = unionFind()
        confClusters[tupleInputs].find(atomNum)  # in case of singletons
        for otherNum, bondType in atomBonds[atomNum]:
          if listInputs in self.posClusterLists[otherNum]:
            confClusters[tupleInputs].union(atomNum, otherNum)
    for tupleInputs, clusters in iter(confClusters.items()):
      for atomLists in clusters.toLists():
        #make a conf for each
        thisConf = self.confNums[-1] + 1
        self.confAtoms[thisConf] = atomLists
        self.confInput[thisConf] = tupleInputs
        self.confNums.append(thisConf)
    #print self.confNums, self.confAtoms, self.confInput
    #that's it, confs have been built

  def _findSets(self):
    '''puts conformations together into sets'''
    self.setToConfs = {}  # maps set numbers to conf lists
    for confNum in self.confNums:
      for tupleInput in self.confInput[confNum]:
        if tupleInput not in self.setToConfs:
          self.setToConfs[tupleInput] = []
        self.setToConfs[tupleInput].append(confNum)
    #print self.setToConfs
    #self.setToConfs contains relevant mapping

  def _assignCoords(self, xyzData):
    '''for each conf (including rigid) find atom positions for each atom'''
    self.outAtoms = 0  # counter to indicate how many there are
    self.outAtomOrigAtom = {}  # maps to original atom numbers from mol2
    self.outAtomInputConf = {}
    self.outAtomConfNum = {}
    self.confNumAtomList = {}
    for confNum in self.confNums:
      self.confNumAtomList[confNum] = []
      for atomNum in self.confAtoms[confNum]:
        self.outAtoms += 1
        globalAtomNum = self.outAtoms
        self.outAtomOrigAtom[globalAtomNum] = atomNum
        self.outAtomInputConf[globalAtomNum] = self.confInput[confNum][0]
        self.outAtomConfNum[globalAtomNum] = confNum
        self.confNumAtomList[confNum].append(globalAtomNum)

  def _identifyClashSetnums(self, clashDecider, mol2data):
    '''for each set decide if it is broken/clashed
    and add it to the self.brokenConfs list if it is. clashDecider is a
    clash.Clash object that figures out what a clash is. mol2data is the
    mol2.Mol2 object that has atom type information and bondedTo method.'''
    self.brokenSets = []
    for aSet in self.setToConfs.keys():
      if clashDecider.decide(mol2data, mol2data.atomXyz[aSet]):
        #means there was a clash
        self.brokenSets.append(aSet)
      #otherwise we do nothing

  def _initClusters(self, clusters):
    '''initializes or reinitializes the clusters of conformations'''
    self.clusters = {}
    self.setNameRemap = {}  # maps old sets to new names
    self.setNameOutOrder = []
    self.setNameFirst = {}
    self.setNameLast = {}
    curSetName = 1
    for clusterIndex, cluster in enumerate(clusters):  # save each cluster
      self.clusters[clusterIndex] = tuple(cluster)
      self.setNameFirst[clusterIndex] = curSetName
      for setName in cluster:
        self.setNameRemap[setName] = curSetName  # map from old to new
        self.setNameOutOrder.append(setName)
        curSetName += 1  # advance counter
      self.setNameLast[clusterIndex] = curSetName - 1  # doing inclusive

  def _findAdditionalMatchSpheres(self, numSpheres=5, cutoff=2.5):
    '''for each cluster, find a couple matching spheres for distant atoms
    that are relatively localized in space.
    data ends up in dict self.clusterSpheres.
    numSpheres is the max# of spheres to add for each cluster. will not always
     find as many as requested.
    cutoff is used as a cutoff to decide
     whether or not to add a sphere for that atom, mean pairwise dist?'''
    atomDists = self.mol2data.distFromAtoms(self.rigidComponent)  # useful
    possibleAtoms = set(self.heavyAtomNums)  # only heavy can be matching
    possibleAtoms.difference_update(self.rigidComponent)  # no need to repeat
    possAtomDist = []  # useful for sorting by distance
    for possibleAtom in possibleAtoms:
      possAtomDist.append((possibleAtom, atomDists[possibleAtom]))
    possAtomDist.sort(key=operator.itemgetter(1), reverse=True)
    #use possAtomDist for each cluster now to find the best candidates
    self.clusterSpheres = {}  # indexed by clusterIndex just like self.clusters
    for clusterIndex in self.clusters.keys():
      cluster = self.clusters[clusterIndex]  # cluster is a tuple of confs
      #print "cluster", cluster #debugging
      self.clusterSpheres[clusterIndex] = []
      for possibleAtom, atomDist in possAtomDist:
        xyzPositions = self.mol2data.getXyzManyConfs(cluster, possibleAtom)
        okayToAdd = False
        if 1 == len(xyzPositions):  # singleton cluster, definitely okay
          okayToAdd = True
        else:
          longDist, meanDist = geometry.longestAndMeanDist(xyzPositions)
          if meanDist <= cutoff:  # passes cutoff
            okayToAdd = True
        if okayToAdd:  # either singleton or passes cutoff
          avgPoint = geometry.getAverage(xyzPositions)
          self.clusterSpheres[clusterIndex].append((possibleAtom, avgPoint))
          #print possibleAtom  # debugging
          if len(self.clusterSpheres[clusterIndex]) == numSpheres:  # done
            break  # out of for loop, no need to go on
      #print self.clusterSpheres[clusterIndex] #debugging

  def _makeClouds(self):
    '''highest level of hierachical ligand sampling, breaks the input
    sets into a few clouds representing gross levels of similar conformations'''
    atomDists = self.mol2data.distFromAtoms(self.rigidComponent)
    #needs switched to divisive bisecting k-means clustering to be fast.
    clusters = self.mol2data.divisiveClustering()
    #printClusterHelper(clusters)  # debug cluster assignments
    self._initClusters(clusters)
    #now that we have clusters, want to find additional matching spheres
    #(with colors even though coloring is bad)
    #data ends up in dict self.clusterSpheres
    self._findAdditionalMatchSpheres()

  def _colorWriter(self, outFile, mol2data):
    '''writes the color table if it was changed from the default'''
    if mol2data.colorConverter.colorInts != \
        mol2data.colorConverter.colorIntsDefault:  # if not default
      colors = mol2data.colorConverter.colorInts.items()
      colors.sort(key=operator.itemgetter(1))
      for colorName, colorKey in colors:
        #outFile.write('T %2d %8s\n' % (colorKey, colorName))
        outFile.write(bytes('T %2d %8s\n' % (colorKey, colorName),'latin1'))

  def _allButSetWriter(
      self, outFile, mol2data, solvdata, setsTotal, clustersTotal=0):
    '''writes the M A B X R and C lines'''
    #now the molecule section, facts about the whole molecule, 5 lines

    outFile.write(bytes(
        'M %16s %9s %3d %3d %6d %6d %6d %6d %6d %6d\n' % (
            mol2data.name[-16:], mol2data.protName[-9:],
            len(mol2data.atomNum), len(mol2data.bondStart),
            self.outAtoms, self.confNums[-1], setsTotal,
            #self.heavyRigidCount, 5, clustersTotal))
            #int(self.heavyRigidCount), int(5), int(clustersTotal)))
            self.heavyRigidCount, 5, clustersTotal),'latin1'))

    #second molecule line, solvation and charge data
    outFile.write(bytes(
        'M %+9.4f %+10.3f %+10.3f %+10.3f %9.3f\n' % (
            solvdata.totalCharge, solvdata.totalPolarSolv,
            solvdata.totalApolarSolv, solvdata.totalSolv,
            solvdata.totalSurface),'latin1'))
    #smiles and long version of name
    outFile.write(bytes('M %-76s\n' % (mol2data.smiles[-76:]),'latin1'))
    outFile.write(bytes('M %-76s\n' % (mol2data.longname[-76:]),'latin1'))
    #best dud energy, computed and put in later. idea is to store the best
    #energy that can be found using the old DOCK/db methods and make sure
    #we aren't totally missing the ball.
    outFile.write(bytes('M %+10.4f\n' % 999.999,'latin1'))
    #atom line, 1 per atom
    for atomNum in range(len(mol2data.atomNum)):
      outFile.write(bytes(
          'A %3d %-4s %-5s %2d %2d %+9.4f %+10.3f %+10.3f %+10.3f %9.3f\n' % (
              mol2data.atomNum[atomNum], mol2data.atomName[atomNum],
              mol2data.atomType[atomNum],
              mol2data.dockNum[atomNum], mol2data.colorNum[atomNum],
              solvdata.charge[atomNum], solvdata.polarSolv[atomNum],
              solvdata.apolarSolv[atomNum], solvdata.solv[atomNum],
              solvdata.surface[atomNum]),'latin1'))
    #now all the bonds.
    for bondNum in range(len(mol2data.bondStart)):
      outFile.write(bytes(
          'B %3d %3d %3d %-2s\n' % (
              mol2data.bondNum[bondNum], mol2data.bondStart[bondNum],
              mol2data.bondEnd[bondNum], mol2data.bondType[bondNum]),'latin1'))
    #now all the coordinates. this section is complex to output since not
    # all atoms*input coordinates are output.
    for xyzNum in range(self.outAtoms):
      xyzNum += 1  # 1-index nonsense
      atomNum = self.outAtomOrigAtom[xyzNum]
      inputConfNum = self.outAtomInputConf[xyzNum]
      confNum = self.outAtomConfNum[xyzNum]
      xyz = self.mol2data.atomXyz[inputConfNum][atomNum]
      #atomnum needs incremented by 1 to make it match up with the input atom#
      outFile.write(bytes(
          'X %9d %3d %6d %+9.4f %+9.4f %+9.4f\n' %
          (xyzNum, atomNum+1, confNum, xyz[0], xyz[1], xyz[2]),'latin1'))
      #amazingly these coordinates are not converted to integers.
    #rigid xyzs, or really just the ligand xyzs to be used for matching
    self.rigidNumSeen = 0
    for rigidNum in self.heavyRigidAtomNums:
      self.rigidNumSeen += 1
      atomColor = mol2data.colorNum[rigidNum]
      xyz = self.mol2data.atomXyz[0][rigidNum]
      outFile.write(bytes(
          'R %6d %2d %+9.4f %+9.4f %+9.4f\n' %
          (self.rigidNumSeen, atomColor, xyz[0], xyz[1], xyz[2]),'latin1'))
    #conformations...
    for confNum in self.confNums:
      coordStart = min(self.confNumAtomList[confNum])
      coordEnd = max(self.confNumAtomList[confNum])
      outFile.write(bytes('C %6d %9d %9d\n' % (confNum, coordStart, coordEnd),'latin1'))

  def _setWriter(self, outFile, mol2data, solvdata):
    '''writes the S lines. no more limit here.'''
    #set conf list S
    if self.clusters is not None:  # if clusters weren't made
      curSets = self.setToConfs.keys()  # this order is fine
      #curSets.sort()
      sorted(curSets)
    else:
      curSets = self.setNameOutOrder
    for outSetNum, curSet in enumerate(curSets):  # all sets
      if self.clusters is not None:  # if clusters weren't made
        outSetNum += 1  # 1 index since it is fortran
      else:
        outSetNum = self.setNameRemap[curSet]
      curConfs = self.setToConfs[curSet]
      totalConfs = len(curConfs)
      if 0 == totalConfs:  # means there are no children, this shouldn't happen
        print ("set", curSet, "has no conformations in it.", curConfs)
        sys.exit(1)
      else:
        totalLines = int(math.ceil(totalConfs / self.coSePerLine))
        lastLineLen = totalConfs % int(self.coSePerLine)
        if 0 == lastLineLen:
          lastLineLen += int(self.coSePerLine)  # correct count when 0
        #the first line that says how many more are coming and has data
        inInput = 0  # mix-n-match
        confEnergy = 999999.999
        outHydro = 3  # mix-n-match
        #this makes the confEnergy a mmff internal energy, ignoring hydroxyls
        #that have been rotated for now.
        confEnergy = mol2data.inputEnergy[curSet] - min(mol2data.inputEnergy)
        outHydro = mol2data.inputHydrogens[curSet]
        brokenSet = 0  # not broken
        if curSet in self.brokenSets:
          brokenSet = 1  # broken
        outFile.write(bytes(
           'S %6d %6d %3d %1d %1d %+11.3f\n' % (
           outSetNum, totalLines, totalConfs, brokenSet, outHydro,
           confEnergy),'latin1'))
        fullLineFormat = 'S %6d %6d %1d'
        for count in range(int(self.coSePerLine)):
          fullLineFormat += ' %6d'
        fullLineFormat += '\n'
        for lineNum in range(totalLines - 1):  # each full line
          outData = [outSetNum, lineNum + 1, self.coSePerLine]
          for count in range(int(self.coSePerLine)):
            outData.append(
                curConfs[lineNum * int(self.coSePerLine) + count])
          outFile.write(bytes(fullLineFormat % tuple(outData),'latin1'))
        #now write last line separately and carefully
        partLineFormat = 'S %6d %6d %1d'
        outData = [outSetNum, totalLines, lastLineLen]
        for count in range(lastLineLen):
          partLineFormat += ' %6d'
          outData.append(
              curConfs[(totalLines - 1) * int(self.coSePerLine) + count])
        partLineFormat += '\n'
        outFile.write(bytes(partLineFormat % tuple(outData),'latin1'))

  def _cloudWriter(self, outFile, mol2data):
    '''write the cloud data'''
    self.cloudNumSeen = 0
    for clusterId in self.clusters.keys():
      outClusId = clusterId + 1
      countSph = len(self.clusterSpheres[clusterId])
      #next line gets around a bug produced when countSph is 0
      maxSphCount = max(self.cloudNumSeen + countSph, self.cloudNumSeen + 1)
      outFile.write(bytes(
          'D %6d %6d %6d %3d %3d %3d\n' % (
              outClusId, self.setNameFirst[clusterId],
              self.setNameLast[clusterId], countSph, self.cloudNumSeen + 1,
              maxSphCount),'latin1'))
      for matchAtom, matchXyz in self.clusterSpheres[clusterId]:
        self.cloudNumSeen += 1  # advance counter
        atomColor = mol2data.colorNum[matchAtom]
        outFile.write(bytes(
            'D %6d %2d %+9.4f %+9.4f %+9.4f\n' % (
                self.cloudNumSeen, atomColor,
                matchXyz[0], matchXyz[1], matchXyz[2]),'latin1'))

  def write(
      self, db2gzFileName, verbose=False, timeit=False,
      limitset=9999999, writeMode='w'):
    '''writes to the new db2 file format. already gzipped.
    writeMode allows append instead of write(over)'''
    try:  # to open the file
      outFile = gzip.GzipFile(db2gzFileName, writeMode)
      try:
        mol2data = self.mol2data
      except AttributeError:
        print ('mol2data missing when output stage encountered.(3)')
        sys.exit(1)
      try:
        solvdata = self.solvdata
      except AttributeError:
        print ('solvdata missing when output stage encountered.(4)')
        sys.exit(1)
      #check if default colors changed, write if they have.
      self._colorWriter(outFile, mol2data)
      self._allButSetWriter(
          outFile, mol2data, solvdata,
          len(self.setToConfs.keys()), len(self.clusters))
      self._setWriter(outFile, mol2data, solvdata)
      if self.clusters is not None:  # if makeclouds was run
        self._cloudWriter(outFile, mol2data)  # this sucks, have to only
           #write clouds for sets that were written. need to rething huge hack
      outFile.write(bytes('E\n','latin1'))  # write the E line here
      outFile.close()
    except IOError:
      print ("error opening output file", db2gzFileName)
      sys.exit(1)
    if verbose:
      print (db2gzFileName + " file written out")

  def writeMol2(
      self, mol2fileName, verbose=False, timeit=False, separateClusters=True):
    '''writes multi-mol2 files instead of db2 files. useful for debugging
    the clustering (or other procedures). each cluster can be written separately
    and will be given a prefix of cluster.00001. etc'''
    if self.clusters is None:
      separateClusters = False  # don't write non-existent clusters
    if separateClusters:
      currentCluster = self.clusters.keys()[0] + 1
      currentPrefix = "cluster." + string.zfill(currentCluster, 5) + "."
      currentName = currentPrefix + mol2fileName
    else:
      currentName = mol2fileName
    try:  # to open the file
      outFile = open(currentName, 'w')
      try:
        mol2data = self.mol2data
      except AttributeError:
        print ('mol2data missing when output stage encountered.(1)')
        sys.exit(1)
      try:
        solvdata = self.solvdata
      except AttributeError:
        print ('solvdata missing when output stage encountered.(2)')
        sys.exit(1)
      if self.clusters is not None:  # if makeclouds was run
        outFile.close()  # close the open and empty file. stupid stupid hack.
        for clusterId in self.clusters.keys():
          currentCluster = clusterId + 1
          currentPrefix = "cluster." + string.zfill(currentCluster, 5) + "."
          currentName = currentPrefix + mol2fileName
          outFile = open(currentName, 'w')
          outNums = []
          for confNumber in range(
              self.setNameFirst[clusterId], self.setNameLast[clusterId] + 1):
            outNum = self.setNameOutOrder[confNumber - 1]  # hate 1-indexing
            outNums.append(outNum)
          self.mol2data.writeMol2File(outFile, outNums)
          if verbose:
            print (currentName + " file written out")
          outFile.close()
      else:
        self.mol2data.writeMol2File(outFile)  # just write them all
      outFile.close()
    except IOError:
      print ("error opening output file", currentName)
      sys.exit(1)
    if verbose:
      print (currentName + " file written out")

#if -1 != string.find(sys.argv[0], "hierarchy.py"):
if -1 != sys.argv[0].find("hierarchy.py"):
  #nothing to do if called from commandline
  pass
