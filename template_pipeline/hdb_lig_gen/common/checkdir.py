#!/usr/bin/env python
'''Check status of dock in this directory.

Michael Mysinger 200705 Created
Michael Mysinger 200909 Added separate cleandir program
Trent Balius 201503 removed FORTRAN STOP check.  
'''

import os
import sys
from optparse import OptionParser

SIZE_LIMIT = 1000000  # 1 Megabyte

# Error types
NOT_SUB = 'not submitted'
SGE_Q = 'in SGE queue'
NO_OUTDOCK = 'error starting DOCK'
RUNNING = 'running now'
DOCK_ERROR = 'error inside DOCK'
SCORED_NONE = 'scored no ligands'
BROKEN = 'incomplete OUTDOCK'
#BROKEN = 'incomplete OUTDOCK, it may still be running.'
DONE_TYPES = [SCORED_NONE]

# based on  http://aspn.activestate.com/ASPN/Cookbook/Python/Recipe/496941
def readlines_backwards(filename, bufsize=8192):
  '''Read a line at a time from the end of the file to the beginning'''
  f = open(filename, 'rb')
  f.seek(0, 2)
  leftover = ''
  pos = f.tell()
  while pos > 0:
    if pos < bufsize:
      bufsize = pos
    f.seek(-bufsize, 1)
    buf = f.read(bufsize) + leftover
    f.seek(-bufsize, 1)
    lines = buf.splitlines(True)
    for oneline in reversed(lines[1:]):
      yield oneline
    leftover = lines[0]
    pos = f.tell()
  yield leftover
  f.close()

def checkdir(indir):
  sub = os.path.join(indir, 'submitted')
  err = os.path.join(indir, 'stderr')
  out = os.path.join(indir, 'OUTDOCK')
  #if not os.path.exists(sub): #i like to clean these up.
  #  return NOT_SUB
  if not os.path.exists(err):
    return SGE_Q
  if not os.path.exists(out):
    return NO_OUTDOCK
  errfile = open(err)
  errdata = errfile.read(SIZE_LIMIT)
  errfile.close()
  # FORTRAN STOP is not printed if GNU compilation
  # or if not submited to the queue.
  if 'FORTRAN STOP' not in errdata:
    if 'at line number' in errdata:
      return DOCK_ERROR
    else:
      return RUNNING
  #if 'at line number' in errdata:
  #    return DOCK_ERROR
  outfile = readlines_backwards(out)
  outdata = next(outfile)
  if not outdata.startswith('elapsed time'):
    return BROKEN
  for i in range(3):
    outdata = next(outfile) + outdata
  if 'no ligands to sort:' in outdata:
    return SCORED_NONE
  return None

def docheckdir(indir='.'):
  '''Test if dock is done.'''
  check = checkdir(indir)
  if check is not None and check not in DONE_TYPES:
    print((indir, check))
    return False
  return True

def main(argv):
  description = "Check status of dock in this diectory."
  usage = "%prog [options]"
  version = "%prog *version 200909* created by Michael Mysinger"
  parser = OptionParser(
      usage=usage, description=description, version=version)
  parser.set_defaults(indir='.')
  parser.add_option(
      "-i", "--indir",
      help="check results inside this directory (default: %default)")
  options, args = parser.parse_args(args=argv[1:])
  if len(args):
    parser.error(
        "program takes no positional arguments.\n" +
        "  Use --help for more information.")
  passed = docheckdir(indir=options.indir)
  return not passed

if __name__ == '__main__':
  sys.exit(main(sys.argv))
