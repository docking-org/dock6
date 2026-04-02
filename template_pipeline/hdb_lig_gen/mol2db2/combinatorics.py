
#Ryan G. Coleman, Brian K. Shoichet 2010
#combinatorics primitives

def allCombinations(inputLists):
  '''takes a list of lists. returns all possible ways of picking one element
  from each list. as a new list.'''
  stack = [[]]
  for inputList in inputLists:
    newStack = []
    for old in stack:
      for inputListItem in inputList:
        oldCopy = old[:]
        oldCopy.append(inputListItem)  # data is added here
        newStack.append(oldCopy)
    stack = newStack
  return stack
