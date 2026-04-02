import os
import sys
import math
#from openeye.oechem import *
import mol2amsol  ## this is a libary Trent Balius and Sudipto Mukherjee wrote. 
import rdkit
from rdkit import Chem

VERBOSE = False

#--------------------------------------------------------------------------------------------------------------------------

def read_ZmatMOPAC(Zmat_file):
    if VERBOSE:
        print ("")
        print ("just entered read_ZmatMOPAC()")
        print ("")

    infile_ZmatMOPAC = open(Zmat_file,'r')
    lines  =  infile_ZmatMOPAC.readlines()

    ZmatMOPAC_lines = {}
    line_key_infile = 0
    line_key_out = 0
    spl = []

    # loop over all the lines in infile_ZmatMOPAC, whose first 3 lines 
    # contain information which can be disregarded
    for line in lines:
         line_key_infile += 1
         #
         if ( line_key_infile < 4): # effect: cut the first three lines off. They must be disregarded.
            pass
         else:
            line_key_out += 1 # after cutting the first 3 lines, the fourth line should be the new first line and so forth.
            #
            # in the first three lines of the Z-matrix, obabel writes 1 too often. This would rise a non-fatal error message in amsol7.1.
            # this error message should be avoided, therefore the first three Z-matrix lines should be slightly modified
            if (line_key_out == 1):
               spl = line.split()
               spl[2] = "0" # instead of 1
               spl[4] = "0" # instead of 1
               spl[6] = "0" # instead of 1
               line = "%-2s %10.6f %2d %11.6f %2d %11.6f %2d %5d %3d %3d\n" % (spl[0],float(spl[1]),int(spl[2]),float(spl[3]),int(spl[4]),float(spl[5]),int(spl[6]),int(spl[7]),int(spl[8]),int(spl[9])) 
            if (line_key_out == 2):
               spl = line.split()
               spl[4] = "0" # instead of 1
               spl[6] = "0" # instead of 1
               line = "%-2s %10.6f %2d %11.6f %2d %11.6f %2d %5d %3d %3d\n" % (spl[0],float(spl[1]),int(spl[2]),float(spl[3]),int(spl[4]),float(spl[5]),int(spl[6]),int(spl[7]),int(spl[8]),int(spl[9])) 
            if (line_key_out == 3):
               spl = line.split()
               spl[6] = "0" # instead of 1
               line = "%-2s %10.6f %2d %11.6f %2d %11.6f %2d %5d %3d %3d\n" % (spl[0],float(spl[1]),int(spl[2]),float(spl[3]),int(spl[4]),float(spl[5]),int(spl[6]),int(spl[7]),int(spl[8]),int(spl[9])) 
               
            ZmatMOPAC_lines[line_key_out] = line
    infile_ZmatMOPAC.close()

    if VERBOSE:
        print ("")
        print ("read_ZmatMOPAC() has finished.")
        print ("")

    return ZmatMOPAC_lines

#-------------------------------------------------------------------------------------------------------------------------

def create_amsol71_inputfile(Path_and_NameZmatMOPACFile,MoleculeName,ZmatMOPAC_Data):

    if VERBOSE:
        print ("")
        print ("just entered the function create_amsol71_inputfile(): ")
        print ("")
        print ("in the case of problems:")
        print ("Make sure that OpenEye OEChem is installed on your system.")
        print ("")

    string_Path_And_NameZmatMOPACFile = Path_and_NameZmatMOPACFile

    # slice off the ending .ZmatMOPAC from string_Path_And_NameZmatMOPACFile by using [0:-10]:

    # open a file for the SM5.42R calculation in water solvent:
    Actual_Amsol71_InputFile_Water = open("%s" % ( string_Path_And_NameZmatMOPACFile[0:-10] + ".in-wat") ,'w')
     
    # open a file for the SM5.42R calculation in hexadecane solvent:
    Actual_Amsol71_InputFile_Hexadecane = open("%s" % ( string_Path_And_NameZmatMOPACFile[0:-10] + ".in-hex") ,'w')

    # the AMSOL7.1 input needs the net charge of the molecule (in its specific protonated state):
    # The net charge will be extracted from temp.mol2 by adding the partial charges in this mol2-file
    # An OpenEye python tool will be used to do this.
    # 
    #mol = OEMol()
    #infile = string_Path_And_NameZmatMOPACFile[0:-10] + ".mol2" 
    #ifs = oemolistream(infile)
 
    #OEReadMolecule(ifs, mol)
    ## OENetCharge():
    ## Determines the net charge on a molecule. If the molecule has specified partial charges, see OEChem's OEHasPartialCharges function,
    ## this function returns the sum of the partial charges rounded to an integer. 
    ## Otherwise this function returns the sum of the formal charges on each atom of the molecule.
    #netcharge = OENetCharge(mol)
    #ifs.close()

    # Trent Balius tried to replaced openeye call with a different function form his and Sudipto's mol2 libary (mol2amsol). 
    #infile = string_Path_And_NameZmatMOPACFile[0:-10] + ".mol2"

    mol2file = string_Path_And_NameZmatMOPACFile[0:-10] + ".mol2"
    smifile = string_Path_And_NameZmatMOPACFile[0:-10] + ".smi"
    mol2file2 = string_Path_And_NameZmatMOPACFile[0:-10] + ".2.mol2"
    ## openbabel will compute mmff94 (gasteiger does not get the right net charge) charges. 
    #comand = "$openbabelPATH/obabel --partialcharge mmff94 -imol2 "+mol2file+" -omol2 -O " + mol2file2
    comand = "$OBABELBASE/bin/obabel --partialcharge mmff94 -imol2 "+mol2file+" -omol2 -O " + mol2file2
    #comand = "$CHEMAXON_PATH/cxcalc -M charge "+mol2file+"-o" + mol2file2
    print(comand)
    os.system(comand)
    ##print(mol2file)
    #mol = mol2amsol.read_Mol2_file(mol2file2)[0]
    #netcharge = int(round(mol2amsol.formal_charge(mol))) 
    #comand = "$openbabelPATH/obabel -imol2 "+mol2file+" -osmi " 
    #of = os.popen(comand)
    #smi = of.readlines()[0].split()[0]
    #print(smi)
    #of.close()
    print(smifile)
    fhsmi = open(smifile,'r')
    smi = fhsmi.readlines()[0].split()[0]
    print(smi)
    fhsmi.close()
    #m = Chem.MolFromSmiles(smi)
    #m = Chem.rdmolfiles.MolFromMol2File(mol2file)
    m2 = Chem.rdmolfiles.MolFromSmiles(smi)
    #netcharge2 = rdkit.Chem.rdmolops.GetFormalCharge(m)
    netcharge = rdkit.Chem.rdmolops.GetFormalCharge(m2)

    #if VERBOSE:
    #    print ("netcharge of molecule in temp.mol2 (sum of partial charges):", netcharge)
    
    #print ("netcharge of molecule in temp.mol2 (sum of partial charges):", netcharge, netcharge2)
    print ("netcharge of molecule in temp.mol2 (sum of partial charges):", netcharge)
    #
    # COMMENT: DOCK's mol2amsol.py could also be used to sum up the partial charges in temp.mol2 file!

    # write the AMSOL7.1 keywords for a SM5.42R point calculation in water to the AMSOL7.1 water input-file:
    Water_Amsol71_SM542R_Keywords = """CHARGE=%s AM1 1SCF TLIMIT=15 GEO-OK SM5.42R\n& SOLVNT=WATER\n""" % netcharge
    Actual_Amsol71_InputFile_Water.write(Water_Amsol71_SM542R_Keywords)

    # write the AMSOL7.1 keywords for a SM5.42R point calculation in hexadecane to the AMSOL7.1 hexadecane input-file:
    Hexadecane_Amsol71_SM542R_Keywords = """CHARGE=%s AM1 1SCF TLIMIT=15 GEO-OK SM5.42R\n& SOLVNT=GENORG IOFR=1.4345 ALPHA=0.00 BETA=0.00 GAMMA=38.93\n& DIELEC=2.06 FACARB=0.00 FEHALO=0.00 DEV\n""" % netcharge
    Actual_Amsol71_InputFile_Hexadecane.write(Hexadecane_Amsol71_SM542R_Keywords)

    # write the name of the currently treated protonated state of the molecule into the AMSOL7.1 file
    # plus the number of atoms in the molecule
    if VERBOSE:
        print ("len(ZmatMOPAC_Data) = number of atoms in molecule : ", len(ZmatMOPAC_Data))
    NumberOfAtomsInMolecule = len(ZmatMOPAC_Data)
    Molecule_Name_NrAtoms = ( "%s %d\n" % (MoleculeName, NumberOfAtomsInMolecule) )
    Actual_Amsol71_InputFile_Water.write(Molecule_Name_NrAtoms)
    Actual_Amsol71_InputFile_Hexadecane.write(Molecule_Name_NrAtoms)

    # write a blank line after the keywords block and the line showing the name of the protonated state of the molecule to the AMSOL7.1 input-files
    blank_line = "\n"
    Actual_Amsol71_InputFile_Water.write(blank_line)
    Actual_Amsol71_InputFile_Hexadecane.write(blank_line)

    # write the lines of the MOPAC Z-matrix to the AMSOL7.1 input-files
    for line_keys in ZmatMOPAC_Data:
        Actual_Amsol71_InputFile_Water.write(ZmatMOPAC_Data[line_keys])
        Actual_Amsol71_InputFile_Hexadecane.write(ZmatMOPAC_Data[line_keys])

    Actual_Amsol71_InputFile_Water.close()
    Actual_Amsol71_InputFile_Hexadecane.close()

    if VERBOSE:
        print ("")
        print ("just finished the function create_amsol71_inputfile(). ")
        print ("")
 
    return

#-------------------------------------------------------------------------------------------------------------------------
def main():
    if VERBOSE:
        print ("")
        print ("just entering main program in make_amsol71_input.py: ")
        print ("")

    if len(sys.argv) != 3: # if no input
        print (" make_amsol71_input.py needs a ZmatMOPAC file as input.")
        print (" The ZmatMOPAC file must have been generated by 'obabel ... -o mopin ...' (MOPAC Internals)")
        print (" Make sure that obabel is installed on your computer !!!")
        return 

    path_and_file_ZmatMOPAC    = sys.argv[1]
    MoleculeName      = sys.argv[2]

    ZmatMOPAC_data = {}
    ZmatMOPAC_data = read_ZmatMOPAC(path_and_file_ZmatMOPAC)

    if VERBOSE:
        print ("")
        print ("calling create_amsol71_inputfile(): ")
        print ("")
    create_amsol71_inputfile(path_and_file_ZmatMOPAC,MoleculeName,ZmatMOPAC_data)

    return
#################################################################################################################
main()
