#!/usr/bin/env python

#Ryan G. Coleman

import string
import sys
import optparse
import mol2
import hydrogens
import time

def mol2hydroxyl(options):
  '''function that reads multi-mol2 input file, rotates or resets hydrogens,
  writes a new multimol2 file back out'''
  if options.timeit:
    timeStart = time.time()
  mol2data = mol2.Mol2(options.mol2file)
  hydrogenRotater = hydrogens.Hydrogens(options.hydrogenfile)
  if options.rotateh or options.reseth:
    hydrogenRotater.findTerminalHydrogens(mol2data)
    if options.verbose:
      print (mol2data.hydrogensToRotate, " terminal hydrogens found.")
    if options.reseth:
      hydrogenRotater.resetHydrogens(mol2data)
    if options.rotateh:
      mol2data = hydrogenRotater.rotateHydrogens(mol2data)
  mol2data.writeMol2(options.mol2outfile)
  if options.timeit:
    timeEnd = time.time()
    print ("time to do everything (total):", timeEnd - timeStart)

if -1 != sys.argv[0].find("mol2hydroxyl.py"):
  description = "Convert multimol2 hydroxyls, rotate, reset, etc"
  usage = "%prog [options]"
  parser = optparse.OptionParser(usage=usage, description=description)
  parser.add_option(
      "-v", "--verbose", action="store_true", dest="verbose", default=False,
      help="lots of debugging output")
  parser.add_option(
      "-a", "--time", action="store_true", dest="timeit", default=False,
      help="timing output")
  parser.add_option(
      "-m", "--mol2", type="string", action="store", dest="mol2file",
      default="db.mol2", help="multimol2 input file, (default: %default)")
  parser.add_option(
      "-o", "--outmol2", type="string", action="store", dest="mol2outfile",
      default="out.mol2", help="mol2 output file, (default: %default)")
  parser.add_option(
      "-y", "--hydrogens", type="string", action="store", dest="hydrogenfile",
      default=None,
      help="terminal hydrogen parameter file, (default: %default)")
  parser.add_option(
      "-r", "--reseth", action="store_true", dest="reseth", default=False,
      help="reset the planar terminal hydrogen positions, (default: %default)")
  parser.add_option(
      "-z", "--norotateh", action="store_false", dest="rotateh", default=True,
      help="rotate the terminal hydrogen positions, (default: %default)")
  #yes you can do reset and rotate in any combination!
  options, args = parser.parse_args()  # default reads from argv[1:]
  if 0 != len(args):
    parser.error("mol2hydroxyl.py takes no positional arguments" +
                 "  Use --help for more information")
  else:
    #print options
    if options.verbose:
      print ("verbose debugging requested.")
      print (options)
    mol2hydroxyl(options)  # main program call
