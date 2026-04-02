#!/usr/bin/env python2.7

#Ryan G. Coleman
#buckets class is used for fast 3d point overlapping
#puts all points into buckets in each dimension
#only compare to reasonably nearby points
#assumes tolerance <<<< 1

from geometry import distL2Squared3
import math

class Bucket3d(object):
  '''buckets class is used for fast	3d point overlapping
  puts all points into buckets in each dimension
  only compare to reasonably nearby points
  assumes tolerance <<<<	1 angstrom'''

  bigBucket = 100  # arbitrary huge bucket size, want to avoid O(n^2)

  def __init__(self, pointList, tolerance):
    '''takes the point list and makes buckets for searching later. O(n)'''
    self.tolerance2 = tolerance ** 2.  # square the tolerance for speed
    self.pointList = pointList
    self.coords = [
        [point[0] for point in pointList],
        [point[1] for point in pointList],
        [point[2] for point in pointList]]
    self.mins = [10000, 10000, 10000]
    self.maxs = [-10000, -10000, -10000]
    for dimension in range(3):
      self.mins[dimension] = int(
          math.floor(min(self.coords[dimension]) - tolerance))
      self.maxs[dimension] = int(
          math.ceil(max(self.coords[dimension]) + tolerance))
    self.buckets = [
        [set() for count in range(1 + self.maxs[0] - self.mins[0])],
        [set() for count in range(1 + self.maxs[1] - self.mins[1])],
        [set() for count in range(1 + self.maxs[2] - self.mins[2])]]
    for dimension in range(3):
      for count, point in enumerate(self.coords[dimension]):
        bucketOne = int(math.floor(point - tolerance)) - self.mins[dimension]
        bucketTwo = int(math.floor(point + tolerance)) - self.mins[dimension]
        for aBucket in range(bucketOne, bucketTwo + 1):
          self.buckets[dimension][aBucket].add(count)
    self.possiblyNearbyPoints = []  # list of sets basically. no order.
    for xCount in range(1 + self.maxs[0] - self.mins[0]):
      for yCount in range(1 + self.maxs[1] - self.mins[1]):
        for zCount in range(1 + self.maxs[2] - self.mins[2]):
          newSet = self.buckets[0][xCount].intersection(
              self.buckets[1][yCount], self.buckets[2][zCount])
          if len(newSet) > 0:
            self.possiblyNearbyPoints.append(list(newSet))

  def getWithinCluster(self, clusters):
    '''souped up for speed version of code. puts nearby points into the
    unionfind data structure 'clusters'. does every possible shortcut i can
    think of for now. super fast now.'''
    #for bucket in self.possiblyNearbyPoints:
    #  print len(bucket),
    #print "bucket lengths"
    for bucket in self.possiblyNearbyPoints:
      #print len(bucket), len(self.pointList)
      indicesLeft = set(range(len(bucket)))
      while len(indicesLeft) > 0:
        oneIndex = indicesLeft.pop()
        oneXyzIndex = bucket[oneIndex]
        if len(bucket) > self.bigBucket:
          thisCluster = clusters.getList(oneXyzIndex)  # O(n), don't do lots
          #print "trying to skip", len(thisCluster), len(bucket)
          if len(thisCluster) >= len(bucket):  # means we should at least quit
            #doing this bucket, nothing left to union
            break
        oneXyz = self.pointList[oneXyzIndex]
        for twoIndex in range(len(bucket)):
          twoXyzIndex = bucket[twoIndex]
          if len(bucket) > self.bigBucket:
            if twoXyzIndex in thisCluster:
              continue  # skip this iteration of the twoIndex for loop
          twoXyz = self.pointList[twoXyzIndex]
          if distL2Squared3(oneXyz, twoXyz) < self.tolerance2:
            clusters.union(oneXyzIndex, twoXyzIndex)
            try:
              indicesLeft.remove(twoIndex)
            except KeyError:
              pass  # really quite okay
        if len(bucket) == len(self.pointList):  # might be able to quit now if
              #all unioned together already after a single pass.
          clusterList = clusters.toLists()
          if len(clusterList) == 1:  # only one cluster means quit now
            return None  # just quit entirely

  def getWithin(self):
    '''returns pairs of points within the tolerance.
    only compare within buckets. slower but doesn't require unionfind
    data structure, kept for testing, etc.'''
    returnPairs = set()
    for bucket in self.possiblyNearbyPoints:
      for oneIndex, oneXyzIndex in enumerate(bucket):
        oneXyz = self.pointList[oneXyzIndex]
        for twoIndex in range(oneIndex + 1, len(bucket)):
          twoXyzIndex = bucket[twoIndex]
          twoXyz = self.pointList[twoXyzIndex]
          if distL2Squared3(oneXyz, twoXyz) < self.tolerance2:
            if twoXyzIndex < oneXyzIndex:
              oneXyzIndex, twoXyzIndex = twoXyzIndex, oneXyzIndex
            returnPairs.add((oneXyzIndex, twoXyzIndex))
    return returnPairs
