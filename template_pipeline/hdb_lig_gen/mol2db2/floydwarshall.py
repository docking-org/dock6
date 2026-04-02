#ryan g. coleman ryangc@mail.med.upenn.edu
# floyd warshall code
#adapted from CLRS of course
#O(n^3) all pairs shortest paths
#just gives distances currently, not actual paths (no Pi matrix)

def makeMatrix(size, infinity=99999999):
  if size > 0:
    retMat = []
    for count in range(size):
      oneRow = [infinity for count in range(size)]
      retMat.append(oneRow)
    return retMat
  else:
    return False

#this version takes a dictionary of neighbors and distances. format is:
# startnode->[[neighbor, dist], [neighbor, dist], [...]]
def floydWarshall(neighbors, infinity=999999999):
  size = len(neighbors)
  oldMat = makeMatrix(size, infinity)
  orderKeys = {}
  #now initialize from the neighbors
  orderedKeys = neighbors.keys()
  #orderedKeys.sort()
  sorted(orderedKeys)
  for order, key in enumerate(orderedKeys):
    orderKeys[key] = order
  for diagonal in range(size):
    oldMat[diagonal][diagonal] = 0
  for key in orderedKeys:
    neighborList = neighbors[key]
    for neigh, dist in neighborList:
      oldMat[orderKeys[key]][orderKeys[neigh]] = dist
      #symmetric case be done later
  newMat = oldMat[:]
  for mac in range(size):
    oldMat = newMat
    for row in range(size):
      for col in range(size):
        newMat[row][col] = min(
            oldMat[row][col], oldMat[row][mac] + oldMat[mac][col])
  return newMat, orderKeys

#testing code
def runTests():
  print ("testing...")
  neighborTest = {}
  neighborTest[0] = [[1, 2]]
  neighborTest[1] = [[0, 2], [2, 1]]
  neighborTest[2] = [[1, 1], [3, 4]]
  neighborTest[3] = [[4, 5], [2, 4]]
  neighborTest[4] = [[3, 5]]
  neighborTest[5] = [[6, 2]]
  neighborTest[6] = [[5, 2]]
  floydWarshall(neighborTest)
