import os
import sys
import math
#from openeye.oechem import *
#import mol2amsol  ## this is a libary Trent Balius and Sudipto Mukherjee wrote. 
import rdkit
from rdkit import Chem

def main():
    mol2file = sys.argv[1]
    smifile  = sys.argv[2]
    name = mol2file
    mol = Chem.rdmolfiles.MolFromMol2File(mol2file)
    smiles = Chem.MolToSmiles(mol)
    print((smiles, name)) 
    fh = open(smifile,'w')
    fh.write('%s %s\n'%(smiles, name))
    #netcharge = rdkit.Chem.rdmolops.GetFormalCharge(mol)
    #print (netcharge)

    return
#################################################################################################################
main()
