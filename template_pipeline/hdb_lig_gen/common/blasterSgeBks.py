#!/usr/bin/env python

#Ryan G. Coleman, Brian K. Shoichet Lab
# modifed by Trent Balius to update to python3, 2020

#import commands
import subprocess
import os

#the code to get the jobid is adapted from
# http://salilab.org/qb3cluster/Chaining_array_jobs_with_SGE

sgeFileHeader = '''#$ -S /bin/csh
#$ -cwd
#$ -p -10
#$ -j yes
#$ -o /dev/null
#$ -e /dev/null
#$ -q all.q

'''

def submit(
    jobScript, headNode='sgehead.bkslab.org', sgeRoot='/opt/sge',
    sgeCell='bks', sgeDir='/opt/sge/bin/lx24-amd64/', qsub='qsub'):
  '''takes an already made script and submits it. tries locally, if that fails
  it ssh's to the head node and submits it. returns the job id for use with
  hold_jid commands in sge'''
  #status, output = commands.getstatusoutput(qsub + ' ' + jobScript)
  status, output = subprocess.getstatusoutput(qsub + ' ' + jobScript)
  if status != 0:  # means the command failed, likely not on a submission node
    #status, output = commands.getstatusoutput(
    status, output = subprocess.getstatusoutput(
        'ssh ' + headNode + ' "setenv SGE_ROOT ' +
        sgeRoot + ' ; setenv SGE_CELL ' + sgeCell + ' ; ' +
        os.path.join(sgeDir, qsub) + ' ' + jobScript + '"')
  print (output)
  tokens = output.split()
  jobid = 0
  for count in range(len(tokens) - 2):
    if 'Your' == tokens[count]:
      jobid = tokens[count + 2].split('.')[0]
  return int(jobid)
