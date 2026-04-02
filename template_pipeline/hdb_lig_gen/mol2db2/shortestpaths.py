#ryan coleman
#adaptation of shortest paths code from my phd for more general case. kind of.
#no objects, just a method, gpl since penn/phd

import priodict

def shortestPaths(nodes, edges, startDist, initialNodes):
  '''simple shortest paths algorithm, using advanced data structure. calculates
  all distances from the starting set of nodes to all other nodes. returns
  nodes to distance mapping. nodes is a list, edges is a dict from nodes to
  other nodes with distance as a tuple, startdist is a float, initialnodes
  is a list.'''
  currentNodes = priodict.priorityDictionary()
  #currentNodes holds data on nodes left to process
  nodeDist = {}
  for initialNode in initialNodes:
    currentNodes[initialNode] = startDist
  while len(currentNodes) > 0:
    currentNode = currentNodes.smallest()
    lastDist = currentNodes.pop(currentNodes.smallest())
    if currentNode not in nodeDist or lastDist < nodeDist[currentNode]:
      #update the dist, add neighbors to heap
      nodeDist[currentNode] = lastDist
      for neighborNode, nbDist in edges[currentNode]:
        newDist = lastDist + nbDist
        if neighborNode not in currentNodes or \
                newDist <= currentNodes[neighborNode]:
          currentNodes[neighborNode] = newDist  # updates prio dict
  return nodeDist
