#!/usr/bin/python
from __future__ import print_function

#################################################################################################################
##
## This program was written by Trent Balius in 
## the Rizzo Research Group at Stony Brook University released in 2012
## 
#################################################################################################################

import math, sys
import os
import mol2

from math import sqrt

#################################################################################################################
# This function reads in a receptor mol2 file and a list of N residues and produces N + 1 grid files.
# N grids for each receptor residue and 1 grid for the remaining residues.
#################################################################################################################

def mkgrid(prefix,list,gridinputfilename):
    list.append("remaining")
    for v in list:
        if (v == 'remaining'): 
            val = str(v)
        elif (int(v) < 10): 
            val = "00"+str(v)
        elif (int(v) < 100):
            val = "0"+str(v)
        elif (int(v) < 1000):
            val = str(v)
        else:
            print (" WARNING: residues are to large. modify if statement.")
            val = str(v)

        mol2prefix = prefix+"_"+str(val)
        print ("create grid for "+mol2prefix+".mol2")
        # DOCKHOME will be replaced during installation; 
        # single quotes to support spaces in pathnames on some OS's.
        gridpath = "'DOCKHOME'/bin/grid"
        os.system("cp "+mol2prefix+".mol2 temp.mol2")
        os.system(gridpath+" -i "+gridinputfilename +" -o temp.out")
        os.system("mv temp.out "+mol2prefix+".out")
        os.system("mv temp.rec.nrg "+mol2prefix+".nrg")
        os.system("mv temp.rec.bmp "+mol2prefix+".bmp")
        os.system("mv temp.rec.grid.mol2 " + mol2prefix+".rec.grid.mol2")
        #sys.exit() 
#receptor_file                  ./temp.mol2
#box_file                       ./box.pdb
#vdw_definition_file            ./vdw.defn
#chemical_definition_file       ./chem.defn
#score_grid_prefix              ./temp.rec
#receptor_out_file              ./temp.rec.grid.mol2
 
    return
    
#################################################################################################################

def on_list(val,list):
   for v in list:
      if val == v:
         return True
   return False

#################################################################################################################

def mkMol2s(file,outputprefix,reslist):
    # pass mol return mol list
    recs = mol2.read_Mol2_file(file)
    
    if (len(recs) != 1):
        print ("rec file does not contain just one mol.")
        sys.exit(0)
    rec = recs[0]

    #Make mol2 files for each residue listed in reslist
    for r in reslist:
        if (r == 'remaining'):
            resid = r
        elif (int(r) < 10):
            resid = "00"+str(r)
        elif (int(r) < 100):
            resid = "0"+str(r)
        elif (int(r) < 1000):
            resid = str(r)
        else:
            print (" WARNING: residues are to large. modify if statement.")
            resid = str(r)


        res_atom_list = []
        res_bond_list = []
        res_residue_list  = {} # dictionary
        for atom in rec.atom_list:
           if (atom.resnum == r):
              res_atom_list.append(atom)
        for bond in rec.bond_list:
           # both bonded atoms are in the specified residue.
           if (rec.atom_list[bond.a1_num - 1 ].resnum == r and rec.atom_list[bond.a2_num - 1].resnum == r ):
               res_bond_list.append(bond);
        #for (residue in rec.residue_list):
        #   if ( residue.resnum == resid ): 
        #       res_residue_list.append(residue)
        res_residue_list[r] = rec.residue_list[r]
        name = "resid" +  str(resid)
        temp = mol2.Mol(name,res_atom_list,res_bond_list,res_residue_list)
        mol2.write_mol2(temp,outputprefix+"_"+str(resid)+".mol2")

    #Make 1 mol2 file containing all residues not in the list
    # (a) Get list of residues from MOL not in reslist.
    allresidues_dic = {};
    allresidues_remainderlist = [];
    for atom in rec.atom_list:
        if atom.resnum not in allresidues_dic:
            if (not on_list(atom.resnum,reslist)):
               allresidues_dic[atom.resnum] = 1
               allresidues_remainderlist.append(atom.resnum)

    # (b) create mol2 of remaining residues 
    res_atom_list = []
    res_bond_list = []
    res_residue_list  = {} # dictionary
    for resid in allresidues_remainderlist:
        for atom in rec.atom_list:
           if (atom.resnum == resid):
              res_atom_list.append(atom)
        res_residue_list[resid] = rec.residue_list[resid]
    for bond in rec.bond_list:
       # both bonded atoms are in the specified residue.
       if (rec.atom_list[bond.a1_num - 1 ].resnum in allresidues_dic and rec.atom_list[bond.a2_num - 1].resnum in allresidues_dic):
           res_bond_list.append(bond);
    name = "remaining_residues"
    temp = mol2.Mol(name,res_atom_list,res_bond_list,res_residue_list)
    mol2.write_mol2(temp,outputprefix+"_remaining.mol2")

    return

#################################################################################################################
def main():
    if (len(sys.argv) < 5):  # sys.argv[0] is the script name.
        print (" This script needs the following:")
        print (" (1) receptor mol2, ")
        print (" (2) output prefix ")
        print (" (3) grid input filename ")
        print (" (4) residue list ")
        print (" # of arguments = " + str(len(sys.argv)-1))
        return

    recfile      = sys.argv[1]
    outputprefix = sys.argv[2]
    gridinfile   = sys.argv[3]

    list = []

    # srb dec 2012
    # Is this suppose to work when the residue list is empty ?
    # I guess no; so I have set the argument test above for at least 1 residue.
    for i in range(4,len(sys.argv)):
        print (i, sys.argv[i])
        list.append(int(sys.argv[i]))

    # cheek that it is monotonic increasing
    oldval = 0
    for i in range(len(list)):
        if not (list[i] > oldval):
           print ("list not monotonic increasing")
           sys.exit(0)
        oldval = list[i]

    res_list = mkMol2s(recfile,outputprefix,list)
    mkgrid(outputprefix,list,gridinfile)
    return 
#################################################################################################################
#################################################################################################################
main()
