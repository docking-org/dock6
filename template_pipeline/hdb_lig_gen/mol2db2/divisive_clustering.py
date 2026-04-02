#!/usr/bin/env python

#Ryan G. Coleman, Brian K. Shoichet Lab
#utility for divisive bisective clustering.
#common use is to cluster positions of ligand atoms
#written generally to cluster lists of points, returns indices into the
#original lists as lists of lists where each sub-list is a cluster

import sys
import pca  # for bisective PCA part of clustering
try:
  import LinearAlgebra  # for error handling. ugh. numeric
except:
  try:
    import numpy.linalg as LinearAlgebra  # error handling for numpy
  except:
    pass
  pass

class SplitZeroError(Exception):
  '''error thrown when the cluster fails to split, i.e. all points project
  to the same point'''

  def __init__(self, indicesSplit):
    '''simple assignment, not used currently'''
    self.indicesSplit = indicesSplit

def getListForIndices(pointListList, indices):
  '''for each index, get the pointList. return newListList'''
  newListList = []
  for index in indices:
    newListList.append(pointListList[index])
  return newListList

def findLongestSubList(clusters):
  '''helper function to find the longest sublist'''
  longestIndex = None
  for index in range(len(clusters)):
    if longestIndex is None:
      longestIndex = index
    elif len(clusters[longestIndex]) < len(clusters[index]):
      longestIndex = index
  return longestIndex

def findOrigSplitIndices(origList, splitIndices):
  '''the orig list is a bunch of indices. splitIndices maps into it. return
  2 lists, one for each splitIndices, post-remapping them onto the original'''
  newSplits = [[] for count in range(len(splitIndices))]
  for splitIndex in range(len(splitIndices)):
    for oneIndex in splitIndices[splitIndex]:
      newSplits[splitIndex].append(origList[oneIndex])
  return newSplits

def divisiveClustering(
    pointListList, numClusters=30, limit=None,
    startClusters=None, verbose=False, overlap=0):
  '''utility for divisive bisective clustering.
  common use is to cluster positions of ligand atoms
  written generally to cluster lists of points, returns indices into the
  original lists as lists of lists where each sub-list is a cluster'''
  if startClusters is None:  # start with some clusters already split
    clusters = [[count for count in range(len(pointListList))]]
  else:
    clusters = startClusters
  while len(clusters) < numClusters:  # until we have enough clusters
    #step 1 is find largest cluster
    biggestClusterIndex = findLongestSubList(clusters)
    biggestCluster = clusters[biggestClusterIndex]
    if 1 == len(biggestCluster):
      break  # no reason to keep dividing, all clusters unique!
    if limit is not None and limit > len(biggestCluster):
      break  # no reason to keep dividing, all clusters small enough
    #get just those clusters points
    clusterToSplit = getListForIndices(pointListList, biggestCluster)
    try:
      splitIndices = pca.findProjectAndSplit(
          clusterToSplit, altSplit=False, overlap=overlap)
      if 0 == len(splitIndices[0]) or 0 == len(splitIndices[1]):
        raise SplitZeroError(splitIndices)
      origSplits = findOrigSplitIndices(biggestCluster, splitIndices)
      del clusters[biggestClusterIndex]  # remove this cluster
      clusters.extend(origSplits)  # add the 2 new clusters
    except LinearAlgebra.LinAlgError:
      if verbose:
        print ("convergence problem during clustering. quitting with ")
        print (len(clusters), " clusters, which should be enough for anybody")
      break  # quit now if we can't converge. means cluster too similar.
    except SplitZeroError:
      if verbose:
        print ("projection/split problem during clustering. quitting with ")
        print (len(clusters), " clusters, which should be enough for anybody")
      break  # quit now if we can't converge. means cluster too similar.
    #print clusters #debugging
    if verbose:
      print ("size of each cluster")
      for cluster in clusters:
        print (len(cluster))
      print (" ")  # more debugging
  #if we get here, all clusters are singletons or numClusters have been reached
  return clusters
