'''General utility functions for python scripts

Michael Mysinger 200707 Created
Michael Mysinger 200803 Change read_splits to use generator
'''

import os
import gzip
import bz2

class MMMError(Exception):
    pass

def read_splits(filename, delim=None, raw_file=False):
    '''Read splits from any whitespace delimited file.'''
    if raw_file:
        f = filename
    else:
        f = open(filename, 'r')
    for x in f:
        yield x.split(delim)
    if not raw_file:
        f.close()

def write_splits(filename, splits, delim='\t', raw_file=False):
    '''Write splits to file delimited by delim.'''
    if raw_file:
        f = filename
    else:
        f = open(filename, 'w')
    for x in splits:
        f.write(delim.join(str(y) for y in x)+'\n')
    if not raw_file:
        f.close()

def flex_open(filename, search=True):
    '''Return compressed or normal file, even altering extensions as needed.'''
    if os.path.exists(filename):
        if os.path.splitext(filename)[1] == ".gz":
            return gzip.GzipFile(filename, 'r')
        elif os.path.splitext(filename)[1] == ".bz2":
            return bz2.BZ2File(filename, 'r')
        else:
            return open(filename, 'r')
    elif search:
        if os.path.exists(filename + '.gz'):
            return gzip.GzipFile(filename + '.gz', 'r')
        elif os.path.exists(filename + '.bz2'):
            return bz2.BZ2File(filename + '.bz2', 'r')
        elif os.path.exists(os.path.splitext(filename)[0]):
            return open(os.path.splitext(filename)[0], 'r')
    raise MMMError("No such file: " + filename)

def read_dirlist(indir):
    dirlist = os.path.join(indir, 'dirlist')
    if not os.path.exists(dirlist):
        raise MMMError("Unable to find dirlist '%s'!" % dirlist)
    f = open(dirlist)
    for i in f:
        dirname = i.strip()
        if not dirname:
            continue
        subdir = os.path.join(indir, os.path.basename(i.strip()))
        if os.path.isdir(dirname):
            subdir = dirname
        elif not os.path.isdir(subdir):
            raise MMMError("Unable to find dirlist entry '%s'!" % dirname)
        yield subdir
    f.close()

def bunch(inlist, n):
    '''Package up list elements in sub-lists n at a time.

    Examples:
    x = [1, 2, 3, 4, 5, 6, 7, 8, 9]
    bunch(x, 1) --> [[1], [2], [3], [4], [5], [6], [7], [8], [9]]
    bunch(x, 2) --> [[1, 2], [3, 4], [5, 6], [7, 8], [9, None]]
    bunch(x, 3) --> [[1, 2, 3], [4, 5, 6], [7, 8, 9]]
    bunch(x, 4) --> [[1, 2, 3, 4], [5, 6, 7, 8], [9, None, None, None]]
    '''
    outlist = [inlist[i:i+n] for i in range(0, len(inlist), n)]
    if outlist:
        outlist[-1].extend([None] * (n - len(outlist[-1])))
    return outlist

def split_sequence(inlist, times):
    '''Splits list by taking alternating elements to form times sublists.

    Examples:
    x = [1, 2, 3, 4, 5, 6, 7, 8, 9]
    split_sequence(x, 2) --> [[1, 3, 5, 7, 9], [2, 4, 6, 8, None]]
    split_sequence(x, 3) --> [[1, 4, 7], [2, 5, 8], [3, 6, 9]]
    split_sequence(x, 4) --> [[1, 5, 9], [2, 6, None],
                                 [3, 7, None], [4, 8, None]]
    '''
    outlist = [group(inlist, times, offset) for offset in range(times)]
    if outlist:
        [outlist[i].extend([None] * (len(outlist[0]) - len(outlist[i]))) for i
            in range(len(outlist))]
    return outlist

def group(inlist, times, offset):
    '''Return offset plus every times element of inlist.

    Examples:
    x = [1, 2, 3, 4, 5, 6, 7, 8, 9]
    group(x, 3, 0] --> [1, 4, 7]
    group(x, 3, 1] --> [2, 5, 8]
    group(x, 2, 0] --> [1, 3, 5, 7, 9]
    '''
    return [inlist[i] for i in range(offset, len(inlist), times)]
