
import math

#a bunch of functions using boltzmann

def occupancyToEnergy(occupancy,m):
  '''uses -kT log occupancy to compute the energy of a residue'''
  return -0.593 * math.log(occupancy) * m
