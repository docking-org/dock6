#ryan g coleman, ryangc@mail.med.upenn.edu
#copyright 2006-7 ryan g coleman, kim sharp crystal.med.upenn.edu
#geometric primitives like distance functions and such

import math
import logging
import numpy
from numpy import linalg
from numpy import subtract
import scipy.spatial.distance
useNumpy = True
useNumeric = False
useScipy = True

#try:  # to use numeric
#    import numpy
#    from numpy import linalg
#    from numpy import subtract
#    import scipy.spatial.distance
#    useNumpy = True
#    useNumeric = False
#    useScipy = True
#except ImportError:  # fallback to numpy if possible
#  useNumpy = False
#  useScipy = False
#  try:
#    import Numeric
#    import Matrix
#    import LinearAlgebra
#    useNumeric = True  # use numeric, if available
#  except ImportError:  # otherwise fallback to hard coded single use code
#    useNumeric = False  # found a simple matrix class in pure python
#    logging.warning("Could not find numpy or numeric")
#    try:
#      import pMatrix
#      #http://aspn.activestate.com/ASPN/Cookbook/Python/Recipe/189971
#    except ImportError:
#      logging.critical("Could not find any matrix library (Numeric, numpy, pMatrix)! Cannot proceed")
#      raise


def distL2(x, y):
  return linalg.norm(numpy.subtract(x, y), ord=2)

def distL2Squared(x, y):
  return scipy.spatial.distance.sqeuclidean(x, y)

distL2Squared3 = distL2Squared


def distL2_old(a, b):
  '''no error checking, very fast, should use everywhere'''
  sum = 0.
  for count in range(len(a)):
    sum += (b[count] - a[count])**2.
  return math.sqrt(sum)  # is this faster than **0.5?

def distL2Squared3_old(a, b):
  '''no error checking, unrolled loop'''
  return (b[0] - a[0])**2. + (b[1] - a[1])**2. + (b[2] - a[2])**2.

def distL2Squared_old(a, b):
  '''no error checking, very fast, should use everywhere, doesn't square root'''
  sum = 0.
  for count in range(len(a)):
    sum += (b[count]-a[count])**2.
  return sum

if not useNumpy and False:
    distL2 = distL2_old
    distL2Squared = distL2Squared_old
    distL2Squared3 = distL2Squared3_old


def dist(a, b, metric='L2'):
  '''a and b should be lists of equal length (any dimension)
  calculates distance needed and returns it (L1,L2,LINF,L2SQUARED).
  these new versions are twice the speed of using list comprehensions.'''
  if metric == 'L2':
    sum = 0.
    for count in range(len(a)):
      sum += (b[count]-a[count])**2.
    return sum**0.5
  elif metric == 'LINF':
    max = 0.
    for count in range(len(a)):
      new = abs(b[count]-a[count])
      if new > max:
        max = new
    return max
  elif metric == 'L2SQUARED':
    sum = 0.
    for count in range(len(a)):
      sum += (b[count]-a[count])**2.
    return sum
  elif metric == 'L1':
    sum = 0.
    for count in range(len(a)):
      sum += abs(b[count]-a[count])
    return sum

def longestAndMeanDist(pts):
  dists = scipy.spatial.distance.pdist(pts, metric='euclidean')
  longestDist = dists.max()
  meanDist = dists.mean()
  return longestDist, meanDist

def longestAndMeanDist_old(pts):
  '''given a list of points, finds the largest distance between any 2. also
  finds mean distance between all pairs. returns both, in that order.'''
  longestDist = 0.
  sumDists, countDists = 0., 0
  for indexOne, ptOne in enumerate(pts):
    for ptTwo in pts[indexOne + 1:]:  # no duplicates, minimal looping
      thisDist = distL2(ptOne, ptTwo)
      longestDist = max(thisDist, longestDist)
      sumDists += thisDist
      countDists += 1
  return longestDist, sumDists/float(countDists)

if not useScipy and False:   # Converting to array takes a long long time. Should be converted once
  longestAndMeanDist = longestAndMeanDist_old

def getAngle(a, b):
  '''helper function for triangle interior, returns angle between two vectors'''
  ab = a[0] * b[0] + a[1] * b[1] + a[2] * b[2]  # all inlined for speed
  aSquared = a[0]**2. + a[1]**2. + a[2]**2.
  bSquared = b[0]**2. + b[1]**2. + b[2]**2.
  #ab = 0. #tons of debugging here
  #aSquared = 0.
  #bSquared = 0.
  #for index in xrange(len(a)):
  #  ab += a[index] * b[index]
  #  aSquared += a[index]**2.
  #  bSquared += b[index]**2.
  return math.acos(
      max(-1., min(1., (ab) / (((aSquared)**0.5)*((bSquared)**0.5)))))

def calcTriAreaList(abc):
  '''uses heron's formula'''
  a, b, c = abc  # unpack
  dists = [distL2(a, b), distL2(b, c), distL2(a, c)]
  s = (dists[0] + dists[1] + dists[2])*0.5
  triArea = (s*(s-dists[0])*(s-dists[1])*(s-dists[2]))**(0.5)
  return triArea

def calcTriArea(a, b, c):  # 3 points in 3d
  '''uses heron's formula'''
  dists = [distL2(a, b), distL2(b, c), distL2(a, c)]
  s = (dists[0] + dists[1] + dists[2])*0.5
  triArea = (s*(s-dists[0])*(s-dists[1])*(s-dists[2]))**(0.5)
  return triArea

def getVector(a, b):
  '''does a-b, returns'''
  return [a[i]-b[i] for i in range(len(a))]

def getNormalVector(a, b):
  '''normal(a-b)'''
  return normalizeVector(getVector(a, b))

def getVector(a, b):
  '''does a-b, returns'''
  return [a[i]-b[i] for i in range(len(a))]

def normalizeVector(vector):
  '''divides each by the total components squared'''
  total = 0.
  for coord in vector:
    total += coord**2.
  total = total**0.5
  newVect = []
  for coord in vector:
    newVect.append(coord/total)
  return newVect

def length(vector):
  '''vector length'''
  total = 0.
  for coord in vector:
    total += coord**2.
  total = total**0.5
  return total

def dot(x, y):
  '''gives dot product of two vectors of any dimension, assumes same length'''
  dot = 0.
  for index in range(len(x)):
    dot += x[index] * y[index]
  return dot

def cross(x, y):
  '''gives cross product of two vectors'''
  return [
      x[1] * y[2] - x[2] * y[1],
      x[2] * y[0] - x[0] * y[2],
      x[0] * y[1] - x[1] * y[0]]

def getDihedralUnited(all):
  '''list of 4 xyzs, gets the dihedral'''
  return getDihedral(all[0], all[1], all[2], all[3])

def getDihedral(a, b, c, d):
  '''4 xyzs, gets the dihedral'''
  cross1 = normalizeVector(
      cross(getNormalVector(a, b), getNormalVector(b, c)))
  cross2 = normalizeVector(
      cross(getNormalVector(b, c), getNormalVector(c, d)))
  try:
    dihedral1 = math.acos(dot(cross1, cross2))
  except ValueError:
    dihedral1 = 0.0  # sometimes the dot ends up a tiny bit above 1.0
  #have to figure out +- direction
  planeD = calculatePlaneD(cross1, b)
  planeFull = (cross1[0], cross1[1], cross1[2], planeD)
  if not checkPlaneSide(planeFull, d):
    dihedral1 = -dihedral1
  return dihedral1

def rotateAboutLine(aIn, dIn, xyz, theta):
  '''rotates the point xyz about the line d-a to an angle of theta radians'''
  #based on http://inside.mines.edu/~gmurray/ArbitraryAxisRotation/
  #  ArbitraryAxisRotation.html
  #first we have to constrain theta to be within -pi to +pi
  while theta < math.pi:
    theta += 2 * math.pi
  while theta > math.pi:
    theta -= 2 * math.pi
  da = getVector(dIn, aIn)  # line through a and d
  #break down and just use the worst notation ever. someone punch me in the face
  a, b, c = aIn  # unpack many things
  d, e, f = dIn
  u, v, w = da
  x, y, z = xyz
  #shortcuts
  uvw = length(da)
  uvw2 = uvw * uvw
  #long stupid equations
  newX = (
      a * (v**2. + w**2.) + u * (- b * v - c * w + u * x + v * y + w * z) +
      (- a * (v**2. + w**2.) + u * (b * v + c * w - v * y - w * z) +
          x * (v**2. + w**2.)) * math.cos(theta) +
      (- c * v + b * w - w * y + v * z) * math.sin(theta) * uvw) / uvw2
  newY = (
      b * (u**2. + w**2.) + v * (- a * u - c * w + u * x + v * y + w * z) +
      (- b * (u**2. + w**2.) + v * (a * u + c * w - u * x - w * z) +
         y * (u**2. + w**2.)) * math.cos(theta) +
      (c * u - a * w + w * x - u * z) * math.sin(theta) * uvw) / uvw2
  newZ = (
      c * (v**2. + u**2.) + w * (- a * u - b * v + u * x + v * y + w * z) +
      (- c * (v**2. + u**2.) + w * (a * u + b * v - u * x - v * y) +
          z * (v**2. + u**2.)) * math.cos(theta) +
      (- b * u + a * v - v * x + u * y) * math.sin(theta) * uvw) / uvw2
  return newX, newY, newZ

def getTriNormalList(united):
  return getTriNormal(united[0], united[1], united[2])

def getTriNormal(a, b, c, firstTime=True):
  '''a, b and c are triange points in clockwise order, returns normal vector
  that points out. returns NORMALIZED vector now. or 0s.'''
  #find a-b and c-b
  #vecAB = normalizeVector(getVector(a, b))
  #vecCB = normalizeVector(getVector(c, b))
  vecAB = getVector(a, b)
  vecCB = getVector(c, b)
  #does the cross product, that's all there is to it
  normal = cross(vecAB, vecCB)
  #only enter this part if all 0 and if first time being called
  if not firstTime:  # has been called recursively.
    return normal  # don't check 0s.don't normalize
  elif firstTime and normal[0] == 0. and normal[1] == 0. and normal[2] == 0.:
    '''this is a big problem. attempt to call after permuting values'''
    newNor = getTriNormal(b, c, a, firstTime=False)  # still maintains clockwise
    if newNor[0] == 0. and newNor[1] == 0. and newNor[2] == 0.:
      lastNo = getTriNormal(c, a, b, firstTime=False)  # again
      #if this is zero we still have to return it
      if lastNo[0] == 0. and lastNo[1] == 0. and lastNo[2] == 0.:
        return lastNo  # 0s knowingly returned
      else:
        return normalizeVector(lastNo)
    else:
      return normalizeVector(newNor)
  else:
    return normalizeVector(normal)

def getAverage(listPoints):
  '''averages any number of 3d points passed in as list'''
  average = [0., 0., 0.]
  for point in listPoints:
    for index in range(len(average)):
      average[index] += point[index]
  for index in range(len(average)):
    average[index] /= len(listPoints)
  return average

def getAverage1(listPoints):
  '''averages any number of 1d points passed in as list'''
  average = 0.
  for point in listPoints:
    average += point
  average /= len(listPoints)
  return average

def getAverageArbitraryDimension(listPoints, dimension=2):
  '''averages any number of nD points passed in as list'''
  average = [0. for count in range(dimension)]
  for point in listPoints:
    for index in range(len(average)):
      average[index] += point[index]
  for index in range(len(average)):
    average[index] /= len(listPoints)
  return average

def planeDistToOrigin(normal):
  '''uses formula from http://mathworld.wolfram.com/Plane.html
  normal is a, b, c, d of plane
  dist = d / ((a^2 + b^2 + c^2) ^ (1/2))'''
  a, b, c, d = normal  # unpack tuple for laziness
  return d / ((a**2. + b**2. + c**2.) ** 0.5)

def fixNormalZeros(vector):
  '''if all 0s, return unchanged, that's fine.
  if 1 or 2 0s, permute a tiny bit so there are no 0s. normalize and return'''
  alpha = 0.0000000000000000001
  if vector[0] == 0. and vector[1] == 0. and vector[2] == 0.:
    return vector  # all zeros
  elif vector[0] == 0. or vector[1] == 0. or vector[2] == 0.:
    newVec = vector[:]  # deep copy, since gets modified
    if vector[0] == 0.:
      newVec[0] += alpha
    if vector[1] == 0.:
      newVec[1] += alpha
    if vector[2] == 0.:
      newVec[2] += alpha
    return normalizeVector(newVec)
  else:
    return vector  # no zeros

def withinTolerance(pointA, pointB, tolerance):
  '''trying to make something fast to check if pointA and pointB are within
  the tolerance of each other.
  exact distance function (l2, l1, linf) not a big deal'''
  if abs(pointA[0] - pointB[0]) < tolerance:
    if abs(pointA[1] - pointB[1]) < tolerance:
      if abs(pointA[2] - pointB[2]) < tolerance:
        return True
  return False

def perturbTriangle(p1, p2, p3):
  '''used to change triangles slightly for intersection checks'''
  p1new = [x+.0000001 for x in p1]
  p2new = [x-.000001 for x in p2]
  p3new = [x+.00001 for x in p3]
  return p1new, p2new, p3new

#p1, p2, p3 are the plane, p4, p5 are the line
#returns the point that is the intersection
#doesn't do uniqueness checks, etc.
#math from Eric W. Weisstein. "Line-Plane Intersection."
#From MathWorld--A Wolfram Web Resource.
#http://mathworld.wolfram.com/Line-PlaneIntersection.html
# t = - |1  1  1  1 |
#       |x1 x2 x3 x4|
#       |y1 y2 y3 y4|
#       |z1 z2 z3 z4|
#       ----------------
#       |1  1  1  0    |
#       |x1 x2 x3 x5-x4|
#       |y1 y2 y3 y5-y4|
#       |z1 y2 z3 z5-z4|
#plug t into:
# x = x4 + (x5-x4)t
# y = y4 + (y5-z4)t
# z = z4 + (y5-z4)t
#uses pMatrix class for now--maybe switch to numericpython if needed
def linePlaneIntersection(p1, p2, p3, p4, p5):
  top = pMatrix.pMatrix(
      [
          [1., 1., 1., 1.],
          [p1[0], p2[0], p3[0], p4[0]], [p1[1], p2[1], p3[1], p4[1]],
          [p1[2], p2[2], p3[2], p4[2]]])
  topDet = top.determinant()
  bottom = pMatrix.pMatrix(
      [
          [1., 1., 1., 0.],
          [p1[0], p2[0], p3[0], p5[0] - p4[0]],
          [p1[1], p2[1], p3[1], p5[1] - p4[1]],
          [p1[2], p2[2], p3[2], p5[2] - p4[2]]])
  botDet = bottom.determinant()
  if topDet == 0.0 or botDet == 0.0:
    return False
  t = -topDet/botDet
  x = p4[0] + (p5[0]-p4[0]) * t
  y = p4[1] + (p5[1]-p4[1]) * t
  z = p4[2] + (p5[2]-p4[2]) * t
  return [x, y, z]

#p1, p2, p3 are the plane, p4, p5 are the line
#returns the point that is the intersection
#doesn't do uniqueness checks, etc.
#math from Eric W. Weisstein. "Line-Plane Intersection."
# From MathWorld--A Wolfram Web Resource.
# http://mathworld.wolfram.com/Line-PlaneIntersection.html
# t = - |1  1  1  1 |
#       |x1 x2 x3 x4|
#       |y1 y2 y3 y4|
#       |z1 z2 z3 z4|
#       ----------------
#       |1  1  1  0    |
#       |x1 x2 x3 x5-x4|
#       |y1 y2 y3 y5-y4|
#       |z1 y2 z3 z5-z4|
#plug t into:
# x = x4 + (x5-x4)t
# y = y4 + (y5-z4)t
# z = z4 + (y5-z4)t
#uses NumericPython for matrix stuff... falls back to pMatrix standalone funct
def linePlaneIntersectionNumeric(p1, p2, p3, p4, p5):
  if not useNumeric:
    return linePlaneIntersection(p1, p2, p3, p4, p5)
  if useNumpy:
    top = [
        [1., 1., 1., 1.],
        [p1[0], p2[0], p3[0], p4[0]], [p1[1], p2[1], p3[1], p4[1]],
        [p1[2], p2[2], p3[2], p4[2]]]
    topDet = numpy.linalg.det(top)
    bottom = [
        [1., 1., 1., 0.], [p1[0], p2[0], p3[0], p5[0]-p4[0]],
        [p1[1], p2[1], p3[1], p5[1]-p4[1]], [p1[2], p2[2], p3[2], p5[2]-p4[2]]]
    botDet = numpy.linalg.det(bottom)
  else:  # actually use numeric
    top = Matrix.Matrix(
        [[1., 1., 1., 1.], [p1[0], p2[0], p3[0], p4[0]], [p1[1], p2[1],
         p3[1], p4[1]], [p1[2], p2[2], p3[2], p4[2]]])
    topDet = LinearAlgebra.determinant(top)
    bottom = Matrix.Matrix(
        [[1., 1., 1., 0.], [p1[0], p2[0], p3[0], p5[0]-p4[0]], [p1[1],
         p2[1], p3[1], p5[1]-p4[1]], [p1[2], p2[2], p3[2], p5[2]-p4[2]]])
    botDet = LinearAlgebra.determinant(bottom)
  if topDet == 0.0 or botDet == 0.0:
    return False
  t = -topDet/botDet
  x = p4[0] + (p5[0]-p4[0]) * t
  y = p4[1] + (p5[1]-p4[1]) * t
  z = p4[2] + (p5[2]-p4[2]) * t
  return [x, y, z]

def intPointInsideTri(p1, p2, p3, intPt):
  '''helper function that checks to see if the intPt is inside
   the triangle p1, p2, p3
  do three checks, make sure intPt is closer to every
  set of 2 vectors than they are to each other'''
  #print "p1, p2, p3, intPt =", p1,",", p2,",", p3,",", intPt
  p2p3ang = getAngle(getVector(p2, p1), getVector(p3, p1))
  if p2p3ang < getAngle(getVector(p2, p1), getVector(intPt, p1)) or \
     p2p3ang < getAngle(getVector(p3, p1), getVector(intPt, p1)):
    return False
  p1p2ang = getAngle(getVector(p1, p3), getVector(p2, p3))
  if p1p2ang < getAngle(getVector(p2, p3), getVector(intPt, p3)) or \
     p1p2ang < getAngle(getVector(p1, p3), getVector(intPt, p3)):
    return False
  p3p1ang = getAngle(getVector(p3, p2), getVector(p1, p2))
  if p3p1ang < getAngle(getVector(p3, p2), getVector(intPt, p2)) or \
     p3p1ang < getAngle(getVector(p1, p2), getVector(intPt, p2)):
    return False
  return True

def intPointInsideTriTuple(triTuple, intPt):
  '''helper function that checks to see if the intPt is inside the
  triangle p1, p2, p3'''
  # the tuple format is ((x), (y), (z), (x-y), (y-x), (y-z), (z-y), (x-z),(z-x))
  #do three checks, make sure intPt is closer to every
  # set of 2 vectors than they are to each other
  inside = True
  #print "triTuple, intPt =", triTuple,",", intPt
  p2p3ang = getAngle(triTuple[4], triTuple[8])
  if p2p3ang < getAngle(triTuple[4], getVector(intPt, triTuple[0])) or \
     p2p3ang < getAngle(triTuple[8], getVector(intPt, triTuple[0])):
    return False
  p1p2ang = getAngle(triTuple[7], triTuple[5])
  if p1p2ang < getAngle(triTuple[7], getVector(intPt, triTuple[2])) or \
     p1p2ang < getAngle(triTuple[5], getVector(intPt, triTuple[2])):
    return False
  p3p1ang = getAngle(triTuple[3], triTuple[6])
  if p3p1ang < getAngle(triTuple[3], getVector(intPt, triTuple[1])) or \
     p3p1ang < getAngle(triTuple[6], getVector(intPt, triTuple[1])):
    return False
  return inside

def getTriNormalList(united):
  return getTriNormal(united[0], united[1], united[2])

def getTriNormal(a, b, c, firstTime=True):
  '''a, b and c are triange points in clockwise order, returns normal vector
  that points out. returns NORMALIZED vector now. or 0s.'''
  #find a-b and c-b
  #vecAB = normalizeVector(getVector(a, b))
  #vecCB = normalizeVector(getVector(c, b))
  vecAB = getVector(a, b)
  vecCB = getVector(c, b)
  #does the cross product, that's all there is to it
  normal = cross(vecAB, vecCB)
  #only enter this part if all 0 and if first time being called
  if not firstTime:  # has been called recursively. don't check 0s.
    return normal  # don't normalize
  elif firstTime and normal[0] == 0. and normal[1] == 0. and normal[2] == 0.:
    '''this is a big problem. attempt to call after permuting values'''
    newNor = getTriNormal(b, c, a, firstTime=False)  # still maintains clockwise
    if newNor[0] == 0. and newNor[1] == 0. and newNor[2] == 0.:
      lastNo = getTriNormal(c, a, b, firstTime=False)  # again
      #if this is zero we still have to return it
      if lastNo[0] == 0. and lastNo[1] == 0. and lastNo[2] == 0.:
        return lastNo  # 0s knowingly returned
      else:
        return normalizeVector(lastNo)
    else:
      return normalizeVector(newNor)
  else:
    return normalizeVector(normal)

def getAverage(listPoints):
  '''averages any number of 3d points passed in as list'''
  average = [0., 0., 0.]
  for point in listPoints:
    for index in range(len(average)):
      average[index] += point[index]
  for index in range(len(average)):
    average[index] /= len(listPoints)
  return average

def getAverageArbitraryDimension(listPoints, dimension=2):
  '''averages any number of nD points passed in as list'''
  average = [0. for count in range(dimension)]
  for point in listPoints:
    for index in range(len(average)):
      average[index] += point[index]
  for index in range(len(average)):
    average[index] /= len(listPoints)
  return average

def findMinsMaxsSpheres(spheres):
  '''goes through all spheres, finds the min and max in each dimension.
  spheres are expected in [x, y, z, r] format'''
  if 0 == len(spheres):
    return False, False  # indicates failure
  mins, maxs = [], []
  for xyz in range(3):
    mins.append(spheres[0][xyz] - spheres[0][3])  # x-radius then y-rad, z-rad
    maxs.append(spheres[0][xyz] + spheres[0][3])  # x+radius then y+rad, z+rad
  for sphere in spheres[1:]:  # already did the first
    for xyz in range(3):
      mins[xyz] = min(mins[xyz], sphere[xyz]-sphere[3])
      maxs[xyz] = max(maxs[xyz], sphere[xyz]+sphere[3])
  return mins, maxs

def lineSphereIntersection(minLine, maxLine, sphere):
  '''line goes from minline to maxline, sphere is x, y, z,radius,
  returns 2 points of intersection, or if failure returns False
  math is from http://en.wikipedia.org/wiki/Ray-sphere_intersection'''
  #move sphere and line so that line starts at 0, 0, 0
  newSphere = []
  for coord in range(3):
    newSphere.append(sphere[coord]-minLine[coord])
  newSphere.append(sphere[3])  # radius
  #convert line to necessary form
  dirLine = []
  for coord in range(3):
    dirLine.append(maxLine[coord]-minLine[coord])
  dirLine = normalizeVector(dirLine)
  partA = 0.
  partB = 0.
  partC = 0.
  for coord in range(3):
    partA += dirLine[coord]*newSphere[coord]  # lxsx + lysx + lzsz
    partB += dirLine[coord]**2.  # lx2 + ly2 + lz2
    partC += newSphere[coord]**2.  # sx2 + sy2 + sz2
  partC -= newSphere[3]**2.  # -sr2
  try:
    oneIntersectionD = (partA + ((partA**2.)-partB*partC)**0.5)/(partB)
    twoIntersectionD = (partA - ((partA**2.)-partB*partC)**0.5)/(partB)
    intersections = [oneIntersectionD, twoIntersectionD]
    if intersections[1] < intersections[0]:
      intersections.reverse()
    #construct output points from original input line
    outputPoints = [[], []]
    for coord in range(3):
      for which in range(2):
        outputPoints[which].append(
            minLine[coord] + dirLine[coord]*intersections[which])
    #print minLine, maxLine, sphere, outputPoints #debugging
    return outputPoints
  except ValueError:
    return False  # didn't work

def countPathTriIntersections(pathPoints, triangle):
  '''checks each line segment against one triangle, counts intersections
  assume pathpoints and triangle have length 3 and are XYZ ordered'''
  intersectionCount = 0
  lastPathPt = pathPoints[0]  # init for loop
  for nextPathPt in pathPoints[1:]:
    triPts0 = triangle[0]
    triPts1 = triangle[1]
    triPts2 = triangle[2]
    posPt, maxIt = False, 5000
    while False == posPt:
      posPt = linePlaneIntersectionNumeric(
          triPts0, triPts1, triPts2, lastPathPt, nextPathPt)
      if False == posPt:
        triPts0, triPts1, triPts2 = perturbTriangle(triPts0, triPts1, triPts2)
        maxIt -= 1
        if maxIt < 0:
          print ("had to perturb points 5000 times", triPts0, triPts1, triPts2, \
              lastPathPt, nextPathPt, "giving up")
          sys.exit(1)
    if posPt is not False:
      if distL2(lastPathPt, nextPathPt) >= distL2(lastPathPt, posPt) and \
          distL2(lastPathPt, nextPathPt) >= distL2(nextPathPt, posPt):
        if intPointInsideTri(triPts0, triPts1, triPts2, posPt):
          # broken  when using large tri?
          intersectionCount += 1
    lastPathPt = nextPathPt  # for next loop
  return intersectionCount

def perturbLine(longAxis, shortAxis1, shortAxis2, startPt, endPt, itersLeft):
    '''makes a slightly different line'''
    #perturb starting line, try again
    newStartPt = [-1., -1., -1.]
    newEndPt = [-1., -1., -1.]
    newStartPt[longAxis] = startPt[longAxis]
    newEndPt[longAxis] = endPt[longAxis]
    if itersLeft % 4 == 3:  # alternate back and forth around line
      newStartPt[shortAxis1] = startPt[shortAxis1] + \
          float(0.0000000001*(5001.-itersLeft))
      newStartPt[shortAxis2] = startPt[shortAxis2] - \
          float(0.000000001*(5001.-itersLeft))
      newEndPt[shortAxis1] = endPt[shortAxis1] + \
          float(0.00000001*(5001.-itersLeft))
      newEndPt[shortAxis2] = endPt[shortAxis2] - \
          float(0.000000001*(5001.-itersLeft))
    elif itersLeft % 4 == 2:  # alternate back and forth around line
      newStartPt[shortAxis1] = startPt[shortAxis1] - \
          float(0.0000001*(5001.-itersLeft))
      newStartPt[shortAxis2] = startPt[shortAxis2] + \
          float(0.000000001*(5001.-itersLeft))
      newEndPt[shortAxis1] = endPt[shortAxis1] + \
          float(0.00000001*(5001.-itersLeft))
      newEndPt[shortAxis2] = endPt[shortAxis2] - \
          float(0.000000001*(5001.-itersLeft))
    elif itersLeft % 4 == 1:  # alternate back and forth around line
      newStartPt[shortAxis1] = startPt[shortAxis1] + \
          float(0.0000000001*(5001.-itersLeft))
      newStartPt[shortAxis2] = startPt[shortAxis2] - \
          float(0.000001*(5001.-itersLeft))
      newEndPt[shortAxis1] = endPt[shortAxis1] - \
          float(0.0000001*(5001.-itersLeft))
      newEndPt[shortAxis2] = endPt[shortAxis2] + \
          float(0.0000000001*(5001.-itersLeft))
    else:
      newStartPt[shortAxis1] = startPt[shortAxis1] - \
          float(0.0000001*(5001.-itersLeft))
      newStartPt[shortAxis2] = startPt[shortAxis2] + \
          float(0.0000001*(5001.-itersLeft))
      newEndPt[shortAxis1] = endPt[shortAxis1] - \
          float(0.000000001*(5001.-itersLeft))
      newEndPt[shortAxis2] = endPt[shortAxis2] + \
          float(0.00000001*(5001.-itersLeft))
    return newStartPt, newEndPt

def getLongestEdge(triList, pointList, direction=-1):
  '''helper function, finds the longest edge in the molecular surface
  direction is 0, 1,2 for the axis to use for projection,
  or -1 to find the euclidean'''
  longestEdge = 0.0
  if -1 == direction:
    for triangle in triList:
      distAB = distL2(
          pointList[triangle[1]-1][1:], pointList[triangle[2]-1][1:])
      distBC = distL2(
          pointList[triangle[2]-1][1:], pointList[triangle[3]-1][1:])
      distCA = distL2(
          pointList[triangle[3]-1][1:], pointList[triangle[1]-1][1:])
      longestEdge = max(distAB, distBC, distCA, longestEdge)
  else:
    pi = [0, 0]
    if 0 == direction:
      pi = [2, 3]  # add 1
    elif 1 == direction:
      pi = [1, 3]  # add 1
    elif 2 == direction:
      pi = [1, 2]  # add 1
    for triangle in triList:
      distAB = distL2(
          [pointList[triangle[1]-1][pi[0]], pointList[triangle[1]-1][pi[1]]],
          [pointList[triangle[2]-1][pi[0]], pointList[triangle[2]-1][pi[1]]])
      distBC = distL2(
          [pointList[triangle[2]-1][pi[0]], pointList[triangle[2]-1][pi[1]]],
          [pointList[triangle[3]-1][pi[0]], pointList[triangle[3]-1][pi[1]]])
      distCA = distL2(
          [pointList[triangle[3]-1][pi[0]], pointList[triangle[3]-1][pi[1]]],
          [pointList[triangle[1]-1][pi[0]], pointList[triangle[1]-1][pi[1]]])
      longestEdge = max(distAB, distBC, distCA, longestEdge)
  return longestEdge

def cacheTriangle(triList, pointList, allowedTris=[-1]):
  '''speed-up function, cache all the various vectors made from a triangle need
  since all triangles get used a couple times, this should be worth it (if you
  have the memory)'''
  #make a vector of tuples
  # the tuple format is ((x), (y), (z), (x-y), (y-x), (y-z), (z-y),
  #                      (x-z), (z-x), (tri#))
  #apparently not [] evaluates to true... so fix that
  cacheDict = {}
  for tri in triList:
    if [-1] == allowedTris or tri[0] in allowedTris:
      x = pointList[tri[1]-1][1:]
      y = pointList[tri[2]-1][1:]
      z = pointList[tri[3]-1][1:]
      xy = getVector(x, y)
      yx = getVector(y, x)
      yz = getVector(y, z)
      zy = getVector(z, y)
      xz = getVector(x, z)
      zx = getVector(z, x)
      tupleRow = (x[0], x[1], x[2]), (y[0], y[1], y[2]), (z[0], z[1], z[2]), \
                 (xy[0], xy[1], xy[2]), (yx[0], yx[1], yx[2]), \
                 (yz[0], yz[1], yz[2]), (zy[0], zy[1], zy[2]), \
                 (xz[0], xz[1], xz[2]), (zx[0], zx[1], zx[2]), \
                 (tri[1], tri[2], tri[3]), (tri[0])
      cacheDict[tri[0]] = tupleRow
  return cacheDict

def calculatePlaneD(normal, pointOnP):
  '''calculates the d of a plane where d = -ax -by -cz where normal = a, b, c
  and point on plane = x, y, z'''
  return - normal[0] * pointOnP[0] - normal[1] * pointOnP[1] - normal[2] * \
      pointOnP[2]

def checkPlaneSide(plane, point):
  '''plane is normal + D (from function calculatePlaneD). sees if point is
  in the direction of normal or not, return boolean'''
  sign = plane[0] * point[0] + plane[1] * point[1] + plane[2] * point[2] + \
      plane[3]
  if sign >= 0:
    return True
  else:
    return False

def planeDistToOrigin(normal):
  '''uses formula from http://mathworld.wolfram.com/Plane.html
  normal is a , b, c, d of plane
  dist = d / ((a^2 + b^2 + c^2) ^ (1 / 2))'''
  a, b, c, d = normal  # unpack tuple for laziness
  return d / ((a**2. + b**2. + c**2.) ** 0.5)

def calculateSphericity(area, volume):
  '''from wikipedia http://en.wikipedia.org/wiki/Sphericity
  Wadell Sphericity, J Geol 1935.
  sphericity = pi^(1/3)(6volume)^(2/3) / area'''
  return ((math.pi**(1./3.))*((6 * volume)**(2. / 3.))) / area
