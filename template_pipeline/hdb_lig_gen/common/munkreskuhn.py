#!/usr/bin/env python

#ryan g. coleman
#implementation of munkres-kuhn algorithm for assignment problem
#based on http://www.public.iastate.edu/~ddoty/HungarianAlgorithm.html
#which is based on:
# Original Reference:  Algorithms for Assignment and Transportation
#Problems, James Munkres, Journal of the Society for Industrial and Applied
#Mathematics Volume 5, Number 1, March, 1957
# Extension to Rectangular
#Arrays Ref:  F. Burgeois and J.-C. Lasalle. An extension of the Munkres
#algorithm for the assignment problem to rectangular matrices.
#Communications of the ACM, 142302-806, 1971.

import sys
import string

def assignAndReturnMatches(cost):
  '''calls assignment, finds the matches, returns a list of them and costs'''
  originalColumns = len(cost[0])
  if originalColumns > len(cost):  # more columns than rows, have to transpose
    return assignAndReturnMatchesTransposed(cost)
  newCost = []  # copy since destroyed
  for line in cost:
    newCost.append(line[:])
  mask = assignment(newCost)
  matches = []  # list of truples
  for row, sublist in enumerate(mask):
    for col, item in enumerate(sublist):
      if 1 == item and col < originalColumns:
        matches.append((row, col, cost[row][col]))
  #print len(matches), len(cost), len(cost[0])
  return matches

def assignAndReturnMatchesTransposed(cost):
  '''TRANSPOSED version for when there are more columns than rows,
  calls assignment, finds the matches, returns a list of them and costs'''
  originalColumns = len(cost)
  newCost = []  # copy since destroyed
  for count in range(len(cost[0])):
    newCost.append([0 for count2 in range(len(cost))])
  for newCol, line in enumerate(cost):
    for newRow, entry in enumerate(line):
      newCost[newRow][newCol] = entry
  #newCost is transposed version of cost!
  #print len(newCost), len(newCost[0])
  #print len(cost), len(cost[0])
  mask = assignment(newCost)
  matches = []  # list of truples
  for row, sublist in enumerate(mask):
    for col, item in enumerate(sublist):
      if 1 == item and col < originalColumns:
        matches.append((col, row, cost[col][row]))
  #print len(matches), len(cost), len(cost[0])
  return matches

def assignment(cost):
  '''sets up run and goes through the six subprocedure steps.
  cost is a list of lists that is used as a 2d matrix.
  returns the set of indices that minimize the cost. should be more rows
  than columns in the cost matrix, or equal.'''
  originalColumns = len(cost[0])
  if originalColumns < len(cost):  # not yet square
    __makeSquare(cost)
  stepNum = 1  # 1 is the first step
  mask = __copyZero(cost)  # 0 = normal, 1 = starred zero, 2 = primed zero
  while 7 != stepNum:  # step 7 is done
    if 1 == stepNum:
      newStepNum = __subtractMinForEachRow(cost)
    elif 2 == stepNum:
      newStepNum = __findZero(cost, mask)
    elif 3 == stepNum:
      rowCover, colCover = __newCovers(cost)
      newStepNum = __findCovers(cost, mask, rowCover, colCover)
    elif 4 == stepNum:
      newStepNum, saveRow, saveCol = __primeOrUncover(
          cost, mask, rowCover, colCover)
    elif 5 == stepNum:
      newStepNum = __constructSeries(
          cost, mask, rowCover, colCover, saveRow, saveCol)
    elif 6 == stepNum:
      newStepNum = __useSmallestValue(cost, mask, rowCover, colCover)
    #print stepNum, newStepNum, cost, mask
    '''
    if stepNum > 3:
      print stepNum, newStepNum
      nonZeroCount,whichOnes = 0, 1
      for rowCount,tempRow in enumerate(mask):
        for colCount,tempEntry in enumerate(tempRow):
          if tempEntry != 0:
            nonZeroCount += 1
            whichOnes *= (rowCount + colCount)
      print nonZeroCount, whichOnes
    '''
    stepNum = newStepNum  # change now for next pass
  return mask

def __makeSquare(matrix):
  '''adds to each row to make columns==rows'''
  rows = len(matrix)
  for row in matrix:
    while len(row) < rows:
      row.append(0)
  #matrix has been modified in place, no need to return

def __copyZero(cost):
  '''makes a copy of the input list of lists in size only, sets all to 0'''
  outMat = []
  for sublist in cost:
    outList = [0 for count in range(len(sublist))]
    outMat.append(outList)
  return outMat

def __subtractMinForEachRow(cost):
  '''finds the minimum in each row (sublist), subtracts it from the sublist.
  cost is modified in-place. always returns 2 as the nextstep. this is step 1'''
  for sublist in cost:
    minval = min(sublist)
    for count in range(len(sublist)):
      sublist[count] -= minval
  return 2  # always go to step 2

def __newCovers(cost):
  '''makes 2 lists of 0s for the dimensions of the fake matrix cost'''
  rowCover = [0 for count in range(len(cost))]
  colCover = [0 for count in range(len(cost[0]))]
  return rowCover, colCover

def __findZero(cost, mask):
  '''finds a zero in a row+column without a starred 0. stars are indicated in
  the mask matrix. uses covers to remember which rows+cols have already been
  checked. modifies cost and mask in place. this is step 2.'''
  rowCover, colCover = __newCovers(cost)
  for row, sublist in enumerate(cost):
    if 0 == rowCover[row]:  # if already marked then we can skip this row
      for col, item in enumerate(sublist):
        if 0. == item and 0 == colCover[col]:  # this can be starred and covered
          mask[row][col] = 1
          rowCover[row] = 1
          colCover[col] = 1
          break  # can skip rest of row since it is covered now
  return 3  # always go to step 3

def __findCovers(cost, mask, rowCover, colCover):
  '''find the covers, i.e. every row and column where there is a starred zero,
  if these cover the row or column of the matrix then you're done, otherwise go
  to step 4. all 4 arguments are modified in place. this is step 3.'''
  for row, sublist in enumerate(cost):
    if 0 == rowCover[row]:  # if already marked then we can skip this row
      for col, item in enumerate(sublist):
        if 1 == mask[row][col]:  # means starred zero
          colCover[col] = 1
          #don't cover rows
  covered = colCover.count(1)
  if covered == len(colCover):  # every column in covered, so we're done
    return 7  # done is step 7
  else:
    return 4  # go to step 4

def __findUncoveredZero(cost, rowCover, colCover):
  '''helper method for primeOrUncover, finds an uncovered zero'''
  for row, sublist in enumerate(cost):
    if 0 == rowCover[row]:  # if already marked then we can skip this row
      for col, item in enumerate(sublist):
        if 0 == cost[row][col] and 0 == colCover[col]:  # means uncovered zero
          return row, col  # done
  return -1, -1  # signal that there is no such uncovered zero

def __isStarInRow(mask, row):
  '''checks the row in question to see if it contains a starred zero'''
  for item in mask[row]:
    if 1 == item:  # starred 0
      return True
  return False  # none found

def __findStarInRow(mask, row):
  '''returns the column of the star in this row, must be star'''
  for col, item in enumerate(mask[row]):
    if 1 == item:
      return col
  return -1  # in case somehow called badly

def __primeOrUncover(cost, mask, rowCover, colCover):
  '''find noncovered zero and prime it. if no starred zero in that row, go to 5.
  otherwise cover the row and uncover the column containing the starred zero.
  continue until no uncovered zeros left. same smallest uncovered value and go
  to 6. this is step 4.'''
  doneYet = False
  while not doneYet:  # infinite loop potential since doneyet is never changed,
                      # however all cases are covered with returns
    uzRow, uzCol = __findUncoveredZero(cost, rowCover, colCover)
    if -1 == uzRow and -1 == uzCol:
      #print "6,-1,-1"
      return 6, -1, -1  # no uncovered zeros, go to step 6, the -1s are fake
    else:  # cover the row and uncover the column containing the starred zero
      mask[uzRow][uzCol] = 2  # prime this zero
      if __isStarInRow(mask, uzRow):
        starCol = __findStarInRow(mask, uzRow)
        rowCover[uzRow] = 1  # cover the row we're in
        colCover[starCol] = 0  # uncover the column of the starred zero
      else:
        return 5, uzRow, uzCol  # go to step 5, save the row+col of the uncover0

def __findStarInColumn(mask, col):
  '''finds a star in this column, returns row, helper for step 5'''
  for row in range(len(mask)):
    if 1 == mask[row][col]:
      return row
  return -1  # failure

def __findPrimeInRow(mask, row):
  '''returns the column of the prime in this row, helper for step 5'''
  for col, item in enumerate(mask[row]):
    if 2 == item:
      return col
  return -1  # failure

def __convertSeries(mask, series):
  '''unstar all stars in series, star all primes, series is list of tuples,
  helper for step 5'''
  for position in series:
    row, col = position
    if 1 == mask[row][col]:  # is a star
      mask[row][col] = 0  # now is unstarred and unprimed
    else:  # must be primed...
      mask[row][col] = 1  # now is starred
  #nothing to return, mask has been changed

def __erasePrimes(mask):
  '''erase all primes (2) in mask, set to 0, helper for step 5'''
  for row, sublist in enumerate(mask):
    for col in range(len(sublist)):
      if 2 == mask[row][col]:  # now is primed
        mask[row][col] = 0  # now is unprimed and unstarred
  #nothing to return, mask has been changed

def __constructSeries(cost, mask, rowCover, colCover, saveRow, saveCol):
  '''construct a series of alternating primed and starred zeros as follows
  z0 = primed zero found in step 4
  z1 = starred zero in the column of z0
  z2 = primed zero in row of z1
  ....
  until primed zero found with no starred zero in its column, then
  unstar each starred zero in series, star each primed zero in series,
  erase primes and uncover everything, go to step 3.
  this is step 5.'''
  zSeries = [(saveRow, saveCol)]  # first set it up
  doneYet = False
  while not doneYet:
    lastColumn = zSeries[-1][1]
    newZrow = __findStarInColumn(mask, lastColumn)
    if -1 != newZrow:  # found it
      zSeries.append((newZrow, lastColumn))
      newZcol = __findPrimeInRow(mask, newZrow)
      zSeries.append((newZrow, newZcol))
    else:  # not found
      doneYet = True  # goes to end after this loop
  __convertSeries(mask, zSeries)
  rowCover, colCover = __newCovers(cost)  # new covers
  __erasePrimes(mask)
  return 3  # go to step 3

def __findSmallestUncovered(cost, rowCover, colCover):
  '''helper for useSmallestValue, finds the smallest uncovered value'''
  minVal = None
  for row, sublist in enumerate(cost):
    if 0 == rowCover[row]:  # if already marked then we can skip this row
      for col, item in enumerate(sublist):
        if 0 == colCover[col]:  # means uncovered column
          if minVal is None or minVal > item:
            #print "item,", item
            minVal = item
  return minVal

def __useSmallestValue(cost, mask, rowCover, colCover):
  '''add smallest uncovered value to every covered row and subtract from every
  column. then go to step 4. this is step 6.'''
  minVal = __findSmallestUncovered(cost, rowCover, colCover)
  #print minVal, rowCover, colCover
  for row, sublist in enumerate(cost):
    for col, item in enumerate(sublist):
      if 1 == rowCover[row]:  # covered row
        cost[row][col] += minVal
      if 0 == colCover[col]:  # means uncovered column
        cost[row][col] -= minVal
  return 4  # go back to step 4

#if called by itself, run tests
if -1 != string.find(sys.argv[0], "munkreskuhn.py"):
  print("assignAndReturnMatches([[1, 2, 3], [4, 5, 1], [3, 1, 10])")
  print((assignAndReturnMatches([[1, 2, 3], [4, 5, 1], [3, 1, 10]])))
  print("  ")
  print(("assignAndReturnMatches(" + \
      "[[1, 2, 3], [4, 5, 1], [3, 1, 10], [50, 20, 30]])"))
  print((assignAndReturnMatches(
      [[1, 2, 3], [4, 5, 1], [3, 1, 10], [50, 20, 30]])))
  print("  ")
  print(("assignAndReturnMatches(" + \
      "[[1, 2, 3], [4, 5, 8], [3, 1, 10], [50, 20, 30]])"))
  print((assignAndReturnMatches(
      [[1, 2, 3], [4, 5, 8], [3, 1, 10], [50, 20, 30]])))
  print("  ")
  print("assignAndReturnMatches([[1, 2, 3], [4, 5, 6], [7, 8, 9]])")
  print((assignAndReturnMatches([[1, 2, 3], [4, 5, 6], [7, 8, 9]])))
  print("  ")
  print(("assignAndReturnMatches(" + \
      "[[1, 2, 3, 4], [4, 5, 6, 7], [7, 8, 9, 10], [11, 12, 13, 15]])"))
  print((assignAndReturnMatches(
      [[1, 2, 3, 4], [4, 5, 6, 7], [7, 8, 9, 10], [11, 12, 13, 15]])))
  print("  ")
  print(("assignAndReturnMatches(" + \
      "[[1, 2, 3, 4], [2, 4, 6, 8], [3, 6, 9, 12], [4, 8, 12, 16]])"))
  print((assignAndReturnMatches(
      [[1, 2, 3, 4], [2, 4, 6, 8], [3, 6, 9, 12], [4, 8, 12, 16]])))
  print("  ")
  print("assignAndReturnMatches([[1, 2, 3], [2, 4, 6], [3, 6, 9]])")
  print((assignAndReturnMatches([[1, 2, 3], [2, 4, 6], [3, 6, 9]])))
