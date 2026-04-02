import sys
from string import *
from math import *

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#calculate cartesian distance between two atoms 
def distance(a, b):
	
	sum = 0
	for i in range(3):
		dist = a[i] - b[i]
		dist = pow(dist, 2)
		sum = dist + sum
	sum = sqrt(sum)
	return sum
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#parse input parameters
cin = sys.argv
r_flag = 'false'
s_flag = 'false'
d_flag = 'false'
for i in range(len(cin)):
	if cin[i] == '-s':
		file1 = cin[i+1]
		s_flag = 'true'
	if cin[i] == '-r':
		file2 = cin[i+1]
		r_flag = 'true'
	if cin[i] == '-d':
		cutoff = atof(cin[i+1])
		c_flag = 'true'
if r_flag == 'false' or s_flag == 'false' or c_flag == 'false':
	print ('ERROR:  missing input parameters!')
	print ('Usage:  python dist_restraints.py -r [receptor.pdb] -s [spheres.sph] -d [distance_cutoff]')


#read in sphere file and collect coordinates
f1 = open(file1, 'r')
sphere_coord = []
while 1:
	next = f1.readline()
	if not next:
		break
	else:
		if next[:4] == 'DOCK':
			continue
		if next[:7] == 'cluster':
			continue
		tmp = split(strip(next))
		sphere_coord.append([atof(tmp[1]), atof(tmp[2]), atof(tmp[3])])
f1.close()

#read in receptor file and collect data
f2 = open(file2, 'r')
pdb = []
while 1:
	next = f2.readline()
	if not next:
		break
	else:
		if next[:4] == 'ATOM' or next[:6] == 'HETATM':
			tmp = split(strip(next))
			pdb.append([atof(tmp[5]), atof(tmp[6]), atof(tmp[7]), tmp[4]])
			
f2.close()

#generate array of flags for receptor atoms
close = []
for i in range(len(pdb)):
	close.append('false')

#determine which receptor atoms are close to spheres
for i in range(len(sphere_coord)):
	for j in range(len(pdb)):
		dist = distance(sphere_coord[i][:3], pdb[j][:3])
		if dist <= cutoff:
			close[j] = 'true'

#collect residues that are close to spheres
res = []
for i in range(len(close)):
	flag = 'lost'
	if close[i] == 'true':
		for j in range(len(res)):
			if pdb[i][3] == res[j]:
				flag = 'found'
		if flag == 'lost':
			res.append(pdb[i][3])

#save resdues in NAB format
out = open('free_residues.txt', 'w')
out.write(':')
for i in range(len(res)-1):
	out.write('%s,' %res[i])
out.write('%s' %res[-1])
out.write('\n')
out.close()
	
		
