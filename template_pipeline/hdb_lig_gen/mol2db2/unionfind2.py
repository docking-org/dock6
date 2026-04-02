#!/usr/bin/env python

#new unionfind code, still using dictionaries but not as weird as other version
#also rewritten from scratch using orig paper + clrs as a guide

class unionFind(object):
  '''allows union, find, tolists, uses path compression and join by rank,
  from CLRS algorithms textbook and previous c implementation by ryan coleman.
  should be much faster than other implementations and w/o object overhead'''

  def __init__(self):
    '''just set up, nothing else done'''
    self._parents = {}
    self._ranks = {}

  def printPar(self):
    '''debugging function'''
    print (self._parents)
    print (self._ranks)

  def check(self, name):
    '''return true if seen already, false otherwise, don't add'''
    if name in self._parents:
      return True
    else:
      return False

  def find(self, name, attachData=False):
    '''looks up name, returns root, does compression along the way'''
    if name not in self._parents:  # hasn't been seen before
      self._parents[name] = name
      self._ranks[name] = 0
      return name
    else:
      path = [name]
      parent = self._parents[name]
      while parent != path[-1]:
        path.append(parent)
        parent = self._parents[parent]
      for item in path[:-1]:  # don't need to do last, already done
        self._parents[item] = parent
        #del self._ranks[item]  # never used again, save space?? XXX??
      return parent

  def union(self, name, other):
    '''unions the 2 sets, inserts them if they haven't been seen before'''
    onePar = self.find(name)
    otherPar = self.find(other)
    if onePar == otherPar:  # already joined
      return onePar
    else:
      if self._ranks[onePar] < self._ranks[otherPar]:  # other is bigger in rank
        self._parents[onePar] = otherPar               # other becomes parent
        return otherPar
      else:                                            # one is bigger in rank
        self._parents[otherPar] = onePar               # one becomes parent
        if self._ranks[onePar] == self._ranks[otherPar]:  # if equal
          self._ranks[onePar] = self._ranks[onePar] + 1   # one's rank goes up
        return onePar

  def different(self, itemA, itemB):
    '''returns true if different, false if same'''
    parA = self.find(itemA)
    parB = self.find(itemB)
    if parA == parB:
      return False
    else:
      return True

  def getList(self, name):
    '''returns the list of all objects in same group as name,
    again this is O(n) since this data is not cached or anything'''
    parent = self.find(name)
    returnList = []
    for item in iter(self._parents.keys()):
      self.find(item)  # makes direct pointer for everything
    for item, par in iter(self._parents.items()):
      if par == parent:
        returnList.append(item)
    return returnList

  def toLists(self):
    '''takes a long time (O(n)), but returns lists of each unioned set'''
    lists = {}
    for item in iter(self._parents.keys()):
      self.find(item)  # makes direct pointer for everything
    for item, par in iter(self._parents.items()):
      try:
        lists[par].append(item)
      except KeyError:
        lists[par] = [item]
    listOfLists = []     # turn into lists of lists
    for aList in iter(lists.values()):
      listOfLists.append(aList)
    return listOfLists

class unionFindAttach(unionFind):
  '''contains attached data carried with each union, unions data together'''

  def __init__(self):
    '''just set up, nothing else done'''
    unionFind.__init__(self)  # call superclass
    self.__attach = {}

  def find(self, name, attachData=False):
    '''looks up name, returns root, does compression along the way'''
    parent = unionFind.find(self, name)  # call superclass
    #print "find", self.__attach, name, attachData
    if parent not in self.__attach:
      self.__attach[parent] = set()
    if attachData and len(attachData) > 0:
      self.__attach[parent].update(attachData)
    #print "findpost", self.__attach, name
    return parent

  def union(self, name, other):
    '''unions the 2 sets, inserts them if they haven't been seen before'''
    onePar = self.find(name)
    otherPar = self.find(other)
    if onePar == otherPar:  # already joined
      return onePar
    else:  # must do union
      if self._ranks[onePar] < self._ranks[otherPar]:  # other is bigger in rank
        self._parents[onePar] = otherPar               # other becomes parent
        self.__attach[otherPar].update(self.__attach[onePar])
        return otherPar
      else:
        self._parents[otherPar] = onePar              # one becomes parent
        if self._ranks[onePar] == self._ranks[otherPar]:  # if equal rank
          self._ranks[onePar] = self._ranks[onePar] + 1   # increment one's rank
        self.__attach[onePar].update(self.__attach[otherPar])
        return onePar

  def getAttached(self, name):
    '''returns the attached set'''
    parent = self.find(name)
    return self.__attach[parent]

  def clearAttached(self, name):
    '''returns the attached set'''
    parent = self.find(name)
    self.__attach[parent] = set()

'''
#commented out this testing code, not even that good of a test
import string,sys #only needed for testing
if -1 != string.find(sys.argv[0], "unionfind2.py"):
  import time
  a = time.time()
  uf = unionFind()
  #uf = unionFindAttach()
  for blah in range(1,5000000):
    #uf.find(blah, set([blah]))
    #print uf.getAttached(blah)
    uf.union(1,blah)
  b = time.time()
  out = uf.toLists()
  c = time.time()
  print b-a, c-b
  #print uf.getAttached(2)
'''
