#!/usr/bin/python
# cp ~jji/bin/md4db-trent.pl 

# Writen by Trent Balius
# this script reads in a file, including files from zinc15.  
# the file will contain the full path to each database file.  
# you can download the file from zinc15 here: 
#    http://zinc15.docking.org/tranches/build/#
# note that the space is very uneven, we will make the chuncks roughly the same size.  

import os
import sys
import math
import shutil

class File_DATA:
    def __init__(self, size, name):
        self.size = size
        self.name = name
    def __cmp__(self, other):
        return cmp(self.size, other.size)
# this defines a compares two LIG_DATA by comparing the two scores
# it is sorted in decinding order.
def bySize(x, y):
    return cmp(x.size, y.size)
def bySizeN(x, y):
    return cmp(-x.size, -y.size)

def is_size_or_count(rtype,size,sizelim,count,countlim):
    #print "in side is_size_or_count"
    #print rtype, size,sizelim,count,countlim
    if (rtype == "size"):
        return (size < sizelim) 
    elif (rtype == "count"):
        return (count < countlim)
    elif (rtype == "both"):
        return ((size < sizelim) and (count < countlim))
    else:
        print "Oh no!"
        return false


if (len(sys.argv) !=  7):
  print len(sys.argv)
  print 'Usage: setup_db2_zinc15_file_number_newstruc.py dir name filename number subnumber type'
  print 'dir: The full path of the directory where docking will be performed, subdirectoriss will be created and contain databases.'
  print 'name: The identifying name of subdirectorys (these are created in the "dir" directory).'
  print 'filename: the file that contains the full paths where dbs are located,\n   this you can download from zinc15.'
  print 'number: number of total diretories to create, this devides the databases into ~ equal docking chunks to be submited.'
  print 'subnumber: number of subdiretories to create, this will make listing files much easer.'
  print 'type: size, count, or both'
  print '      size will make dirs equal chunks in terms of size,'
  print '      "count" option will put equal numbers of database files in the number directory spesified, '
  print '      the option "both" will try and satisfy both size and count options.  '
  sys.exit(-1)

root = "/nfs/db/dockable/zinc15/byproperty";
# this is "Shards"
workdir    = sys.argv[1] # directory to make
name   = sys.argv[2] # 
infile   = sys.argv[3]
tot_number = int(sys.argv[4]) # number of dir, dock chunks to be submited. 
sub_number = int(sys.argv[5]) # number of dir, dock chunks to be submited. 
runtype   = sys.argv[6]

out_number = math.ceil(tot_number / sub_number )

print "workdir = %s\nname = %s\ninfile = %s\ntot_number = %d\nsub_number = %d\nruntype = %s\nout_number = %d\n"%(workdir,name,infile,tot_number,sub_number,runtype,out_number)

# consider making it smarter to determin how meny nested directories are needed. 
if (out_number > sub_number):
    print "tot = %d, subdir = %d, outerdir = %d\n"%(tot_number,sub_number,out_number) 
    print "error. out_number > sub_number there are too meny top directorys"
    exit


fastflag = False 
if (runtype == "count"):
   fastflag = True 

if (runtype == "size"):
   print "type = size"
   print "number (" +str(tot_number)  +  ") is ignored. "

if not (runtype == "size" or runtype == "count" or runtype == "both"):
   print "type = size, count, or both"
   exit()

indock      = "INDOCK"
dirlist     = "dirlist"
gridssphere = "dockfiles"

infh = open(infile,'r')
fileslines = infh.readlines()
infh.close()

#if not (os.path.exists(workdir+"/"+indock+".ori")):
#  os.system("cp "+workdir+"/"+indock+" "+workdir+"/"+indock+".ori")
#else:
#  print workdir+"/"+indock+".ori exists"
#  exit()
#
#os.system('sed -i "s/\.\.\//..\/..\//g" '+ workdir+"/"+indock)

# populate a data structure to sort on size.
# this is one way to do this. 
filedata = []
for line in fileslines:
    #print line
    line = line.strip('\n')
    if not (os.path.exists(line)):
        print line + " doesn't exist."
        splitline =line.split('/')
        line = '' 
        N = len(splitline)
        for i in range(N):
            if i == N-1:
              line = line+'/db2'
            line = line+'/'+splitline[i]
        if not (os.path.exists(line)):
            print line + " doesn't exist."
            continue 
        else: 
            print "the file is "+ line

    if (fastflag):
        #print line
        temp_fd = File_DATA(0, line)
    else: # if not fast then sort by file size. 
        handel = os.popen("ls -lL "+line)
                         #handel = os.popen("ls -1 "+root+"/"+first[i1]+second[i2]+"/"+third+fourth[i4]+fifth[i5]+sixth[i6]+"/db2/*.db2.gz")
        lines = handel.readlines()
        if len(lines) > 1:
            print lines
            exit()
        line = lines[0]
        splitline = line.split()
        
        #print splitline[0], splitline[1], splitline[2], splitline[3], splitline[4]
        temp_fd = File_DATA(int(splitline[4]), splitline[8])

    filedata.append(temp_fd)

numdb  = len(filedata) # number of existing databases, not just all that are in the database index file (some might not exist). 
#numdb  = len(fileslines) # line in file, number of database file to be docked 

print "lines in file = %d;  existing databases=%d; specified dir = %d;\n"%(len(fileslines), len(filedata), int(tot_number)) 

if (numdb < int(tot_number)):
   print "existing dbs in file (%d) < number of directories (%d)."%(numdb,int(tot_number))
   tot_number = numdb

numberchunk = int(numdb) / int(tot_number)
remainder   = int(numdb) % int(tot_number) # this is the number of remaining points, we will distrbute them to the frist set of chunks 

print numberchunk, remainder

# sort by file size, bigest file is at the end of the list.
#filedata.sort(bySize)
# sort by file size, bigest file is at the begining of the list.
if not (fastflag): 
   filedata.sort(bySizeN)
#handel.closed

total = 0 # loop over the files from smallest to largest. 
          # when the size exceds the max make new file.
count = 0

#maxfilesize = filedata[-1].size # ge the value at the end of the list

maxfilesize = filedata[0].size # get the value at the start of the list
count_db = 0
countdir = 0
countlim = numberchunk

outerdircount = 1

#filename = '%s%04d'%(name,count)
filename = 'split_database_index'
#dirname   = '%s/%s%04d'%(dir,name,count)
tdirname   = '%s/t%05d'%(workdir,outerdircount) # this is the top dir name
os.system("mkdir "+tdirname+"/") # make a symbolic link to the dockfiles here. 
#os.system("ln -s "+workdir+"/dockfiles/ "+tdirname) # make a symbolic link to the dockfiles here. 
os.system("ln -s ../dockfiles/ "+tdirname) # make a symbolic link to the dockfiles here. 
dirname   = '%s/t%05d/%s%04d'%(workdir,outerdircount,name,count)

os.system('rm -rf '+dirname+';mkdir -p '+dirname)
#os.system('mkdir -p '+dirname)

if not (os.path.exists(workdir)):
    print workdir + " does not exists.  "
    exit()
elif not (os.path.exists(workdir+"/"+indock)):
    print workdir+"/"+indock + " does not exists.  " 
    exit()
elif not (os.path.exists(workdir+"/"+gridssphere)):
    print workdir+"/"+gridssphere
    exit()
filehdirlist = open(workdir+'/'+dirlist,'w')
filehdirlist.write(dirname+'\n')
fileh = open(dirname+'/'+filename,'w')
#shutil.copy(dir+"/"+indock, dirname)
os.system("cp "+workdir+"/"+indock+" "+dirname)

lfd = len(filedata)
tenpercent = int(round(0.10 * lfd))-1

print lfd, tenpercent

print "----------------------------------------"

per = 10
count_tp = 0
subdircount = 1
for data in filedata:
    #if total maxfilesize:
    if (tenpercent == count_tp or count_db == lfd):
        print '%3d'%per,
        sys.stdout.flush()
        per = per+10 
        count_tp = 0
    else:
        count_tp = count_tp+1

    if is_size_or_count(runtype, total, maxfilesize, count_db, countlim):
        #print data.size, data.name
        fileh.write(data.name+'\n')
    else: 
        #print total
        #print "start new file" 
        #print data.size, data.name
        if (subdircount >= sub_number):
           subdircount = 0
           outerdircount = outerdircount+1
           tdirname   = '%s/t%05d'%(workdir,outerdircount) # this is the top dir name
           os.system("mkdir "+tdirname+"/") # make a symbolic link to the dockfiles here. 
           #os.system("ln -s "+workdir+"/dockfiles/ "+tdirname+"/") # make a symbolic link to the dockfiles here. 
           os.system("ln -s ../dockfiles/ "+tdirname+"/") # make a symbolic link to the dockfiles here. 
        total = 0
        count_db = 0
        fileh.close()
        count = count+1
        #filename = 'temp%04d'%count
        #dirname   = '%s/%s%04d'%(workdir,name,count)
        dirname   = '%s/t%05d/%s%04d'%(workdir,outerdircount,name,count)
        os.system('rm -rf '+dirname+';mkdir -p '+dirname)
        os.system("cp "+workdir+"/"+indock+" "+dirname)
        filehdirlist.write(dirname+'\n')
        fileh = open(dirname+'/'+filename,'w')
        fileh.write(data.name+'\n')
        countdir = countdir + 1
        if countdir > (tot_number-remainder-1): # distribut to chunks
           countlim = numberchunk + 1
           #print "remainder is add to this chunk" 
           #print "countlim == " + str(countlim)
        else: 
           countlim = numberchunk
        subdircount = subdircount+1
    total = total + data.size
    count_db = count_db + 1



#fileh.write(data.name+'\n')
fileh.close()
filehdirlist.close()

print ""
print "----------------------------------------"
print "                                              "
print "                ^^^^^^^^^^^^^^                "
print "               ^^^^^^^^^^^^^^^^               "
print "              ^^^^^^^^^^^^^^^^^^              "
print "             |     ---     ---  \             "
print "  nn   nn    |     <o>   \  <o>  |   LIVE    "
print "   \\\ ///   {            \       \   LONG    "
print "    \\V//|\   |         <..>      |   AND     "
print "   |    | |  |   \___________/|  |   PROSPER "
print "   \      |  \   \  \ |_|_| \ / /        "
print "    \  _ /    \   \__\_|_/__// /          "
print "    /   /    /  \__          _/          "
print "   /   /    /     \_________/            "
print "  /   /    /              /             "
print "----------------------------------------"
