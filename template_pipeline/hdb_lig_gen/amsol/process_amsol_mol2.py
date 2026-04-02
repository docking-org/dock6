import mol2amsol  ## this is a libary Trent Balius and Sudipto Mukherjee wrote. 
import math
import sys
import os
import os.path
import gzip
from math import sqrt

#################################################################################################################
# written by Trent E Balius, B.K. Shoichet lab, Nov. 2013
# and modified and commented with respect to the changes in AMSOL7.1 output-files
# by Thomas B. Adler, B.K. Shoichet lab
#
# this script reads in AMSOL7.1 output-files (for water (.o-wat) and hexadecane (.o-hex))
# and processes them.
# 
# process_amsol_mol2.py is a replacement for the following scripts:  
#   $DOCK_BASE/etc/SubstrSolv2.pl
#   $DOCK_BASE/etc/3Step.csh2 temp.o-wat
#   $DOCK_BASE/etc/3Step.csh2 temp.o-hex
#   $DOCK_BASE/etc/SubstrSolv2.pl temp.o-hex.d temp.o-wat.d temp.solv temp.err
#   $DOCK_BASE/etc/UpdatChrg.pl temp.mol2 temp.solv temp.nmol2 temp.err2
#
# !!! consider adding in: running amsol itself in future.  
#
#################################################################################################################

VERBOSE = False

def is_int(a):
    """Returns true if a can be an integer"""
    try:
        int (a)
        return True
    except:
        return False


#################################################################################################################
#
def process_amsol_file(file,outputprefix,watorhex):
    # reads in data amsol output.
    if VERBOSE:
        print ("")
        print ("**** starting the function process_amsol_file() ****")
        print ("")

    if not (os.path.exists(file)):
        print (file + "does not exist. \n\n Exiting script . . .")
        exit()

    ## read in both gzip and uncompresed file.
    splitfile = file.split('.')
    N = len(splitfile)-1
    if VERBOSE:
        print (splitfile[N])
    if (splitfile[N] == 'gz'):
       file1 = gzip.open(file, 'rb')
    else:
       file1 = open(file,'r')

    ## open up output file
    outputfilename = outputprefix +watorhex +str(".log") 
    file2 = open(outputfilename,'w')
    lines  =  file1.readlines()

    name = ''
    numatoms = 0
    alist = []
    total_line = []

    ## loops over the file
    for line in lines:
         linesplit = line.split() #split on white space
         #  
         # identifying and writing certain lines of the AMSOL7.1 output file into outputfilename.log :
         #
         # infos originally coming from the AMSOL7.1 input, which are reprinted in the AMSOL7.1 output-file:
         # i.e., molecule's name and number of atoms in this molecule:
         # i.e., ZINC00000000.0.1_1 37
         if (len(linesplit) == 2 and is_int(linesplit[1]) and name == ''): # get the name and atom cout.
             file2.write(line)
             name     = linesplit[0]
             numatoms = int(linesplit[1])

         # evaluating the large table close to the end of the AMSOL7.1 output-file :

         # amsol-mod4: table had just 8 columns:
         # if (len(linesplit) == 8 and is_int(linesplit[0]) ):
         # AMSOL7.1: table has an additional column "Sigma k cal/(Ang**2)"
         # AMSOL7.1 table: table has 9 columns
         #if (len(linesplit) == 9 and is_int(linesplit[0]) and ( linesplit[4] is not "*")):  ## this gets that per-atom break-down of solvation calculation
         if (len(linesplit) == 9 and is_int(linesplit[0]) and ( linesplit[4] != "*")):  ## this gets that per-atom break-down of solvation calculation
             file2.write(line)
             alist.append(linesplit)

             # alist is a list running over the atoms i, where each atom i is associated with a list of data:
             # alist contains the following data :::
             # for each atom i, a list of 9 values exists:
             # 
             # alist[i][0] Atom number    (1 up to number of atoms in molecule)
             # alist[i][1] Chem. symbol
             # alist[i][2] CM2 chg.       (CM2 partial atomic charge)
             # alist[i][3] G_P (kcal/mol) (atomic polar contribution to solvation free energy)
             # alist[i][4] Area (Ang**2)  (surface area)
             # alist[i][5] Sigma (kcal/Ang**2) (sigma, cf. "Modeling free Energies of Solvation and Transfer", Computational Thermochemistry, Chapter 15,
             #                                  D.J. Giesen, D.G. Truhlar, 1998)
             # alist[i][6] SS G_CDS (kcal/mol) (atomic apolar contribution to solvation free energy)
             # alist[i][7] Subtotal (kcal/mol) (atomic solvation free energy (polar+apolar))
             # alist[i][8] M value


         ## Extracting the first "Total: ..." line from the amsol-mod4 output-file :
         #elif (len(linesplit) == 6 and linesplit[0] == "Total:" ):
         #
         # Extracting the first "Total: ..." line from the AMSOL7.1 output-file :
         # 6 is not changed, the new column in the large AMSOL7.1 output table is not considered in the first "Total: ..." line 
         elif (len(linesplit) == 6 and linesplit[0] == "Total:" ):
             total_line = linesplit
             # total_line is a list containing now:
             #
             # total_line[0] = "Total:"
             # total_line[1] = CM2 chg.
             # total_line[2] = G_P, i.e. polar contribution to solvation free energy or also called Polarization free energy (in kcal/mol)
             # total_line[3] = Area (Ang**2) (in amsol-mod4: Area (CD) (Ang**2))
             # total_line[4] = SS G_CDS in kcal/mol (in amsol-mod4: G-CD)
             # total_line[5] = Subtotal in kcal/mol, i.e. sum of all atomic polar (G_P[i]) and apolar (SS G_CDS[i]) contributions to the solvation free energy
             #            Subtotal = solvation free energy
             #            Subtotal = "(5)  G-P-CDS(sol) = G-P(sol) + G-CDS(sol) = (2) + (4)             -XX.XXX kcal" 
             #            With respect to Subtotal, be aware of :::
             #            **** NOTA BENE ****
             #            This is the net solvation energy for this exact molecular structure
             #            (nuclear and electronic)!  The standard-state solvation energy should be
             #            obtained as the difference between the heat of formation plus delta-G solvation
             #            for the relaxed solvated system and that for the relaxed gas-phase system.
            
             # in the entire process_amsol_mol2.py machinery :::
             #
             # total_line --> tot_wat --> totwat
             # total_line --> tot_hex --> tothex
             
             file2.write(line)
    file2.close()

    if VERBOSE:
        print ("")
        print (" **** The function process_amsol_file() was finished. ****")
        print ("")

    return alist,total_line,name,numatoms

#################################################################################################################
#################################################################################################################

def diff_amsol71_files(atom_listwat,totwat,atom_listhex,tothex,name,numatom,outputprefix):

    if VERBOSE:
        print ("")
        print ("**** starting the function diff_amsol71_files() ****")
        print ("")
    ## this function will compare hex with wat and write an *.solv output file:
    ## The difference "water minus hexadecane" will be calculated.
    ##
    ## The content of the *.solv file will have the following format:
    ##
    ## ZINC000000  52  0.0   -16.55   436.86     7.27    -9.28
    ## -0.1645   -0.19  19.39    0.56    0.37
    ## -0.6804    0.77   5.64   -0.98   -0.21
    ## -0.3203    0.27   6.50   -0.84   -0.57
    ## -0.1644   -0.92  13.19    0.70   -0.22 
    ##    .        .      .       .       .
    ##    .        .      .       .       .
    ##    .        .      .       .       .
    ##    .        .      .       .       .
    ## 
    ## The first line from the above-mentioned table contains the following information :
    ## ZINC000000  52  0.0   -16.55   436.86     7.27    -9.28 means :
    ##
    ## column 1 molecule name
    ## column 2 number of atoms
    ## column 3 formal charge, i.e., net charge of the entire molecule
    ## column 4 difference of total polarization free energies for water and hex
    ## column 5 Area in (Ang**2) just from hexadecane AMSOL7.1 output-file
    ## column 6 difference of SS G-CDS values from water and hexadecane AMSOL7.1 outputs
    ## column 7 sum of column 4 and column 6: sum of the differences of the total polar and apolar atomic contributions to the solvation free energies for water and hexadecane 
    ##
    ## the following lines are the per-atom break-down:
    ## cf.: 
    ## -0.1645   -0.19  19.39    0.56    0.37 et cetera linea:
    ##
    ## column 1 is atomic partial charge just from the hexadecane (!!!) output file (hexadecane should simulate a protein-like environment: The ligand will be docked
    ##          into a protein pocket. The later use of partial charges obtained for a hexadecane environment in the electrostatic docking term is reasonable.)
    ## column 2 is the difference of the atomic polarization free energies obtained in water and hexadecane(kcal/mol): (wat - hex) 
    ## column 3 is Area (Ang**2) just from hex (is the same for water and hex. hence the source-file (whether hex or water output-files) does not really matter.)
    ## column 4 is SS G-CDS (kcal/mol) (((in amsol-mod 4 called "G-CD")))
    ## column 5 is Subtotal: (kcal/mol) (((in amsol-mod4 output called "Total Solv. free energy")))

    N = len(atom_listwat)
    if len(atom_listwat) != len(atom_listhex):
        print ("Error: len(atom_listwat) != len(atom_listhex):" + str(len(atom_listwat))+" != "+str(len(atom_listhex)))
        sys.exit()

    if N != numatom:
        print ("\n".join(' '.join(l) for l in atom_listwat))
        print ("Error: len(atom_listwat) != numatom: " + str(N) +" != "+str(numatom))
        sys.exit()

    fileline = ''## generate string output to write to a file
    if VERBOSE:
        print ("atom_listwat[i][j] + -- + atom_listhex[i][j] +;   ::::")
        for i in range(N):
            # amsol-mod4 for j in range(0,8):
            for j in range(0,9):
               print (atom_listwat[i][j] +" -- "+ atom_listhex[i][j] +';  ',)
            print ("\n",)

    # initialize the lists/arrays:
    Chghex  = []; Polhex  = []; SAA     = []; 
    Apolhex = []; Polwat  = [];
    Apolwat = [];
    Sigmahex = []
    diff_Pol = []; diff_Apol = []; diff_AtomicSolv = [];

    for i in range(N):
        Chghex.append(0.0)   # list of partial atomic charges from AMSOL7.1 SM5.42R solvation calculation in hexadecane solvent
        Polhex.append(0.0)   # list of atomic polar contributions to solvation free energy: G_P obtained for hexadecane solvent
        SAA.append(0.0)      # list of atomic surface area contributions taken from hexadecane output
        Sigmahex.append(0.0) # list of sigma coefficients 
        Apolhex.append(0.0)  # list of atomic apolar contributions to solvation free energy: SS G_CDS obtained for hexadecane solvent

        Polwat.append(0.0)   # list of atomic polar contributions to solvation free energy: G_P obtained for water solvent
        Apolwat.append(0.0)  # list of atomic apolar contributions to solvation free energy: SS G_CDS obtained for water solvent

        diff_Pol.append(0.0)     # list of differences of the atomic polar contributions to solvation free energy in water and hexadecane: wat - hex
        diff_Apol.append(0.0)    # list of differences of the atomic polar contributions to solvation free energy in water and hexadecane: wat - hex 
        diff_AtomicSolv.append(0.0) # list of differences of the atomic solvation free energies (polar + apolar contributions) in water and hexadecane: wat - hex

    # Initializing the sums (over the atomic contributions):

    Chghexsum = 0.0  # sum of the partial charges (CM2: Truhlar's charge model 2 charges) of each atom in hexadecane (!) solvent
    Polhexsum = 0.0  # sum of the atomic polar contributions to solvation free enthalpy: Polarization Free Energy (G_P) 
    SAAsum = 0.0     # sum of the atomic surface area contibutions in Angstrom
    sum_Apolhex = 0.0    # sum of the atomic apolar contributions to solvation free enthalpy: SS G_CDS (in amsol-mod4: just called G_CD, but contained S)
                   # CDS means: cavity-dispersion-solvent structure (reordering)
                   # sum_Apolhex is the same for water or hexadecane solvent. Therefore, just the hexadecane case has been considered here.
 
    FreeEnergyHex_sum = 0.0       # sum of the atomic contributions to the total solvation free energy obtained in hexadecane: In AMSOL7.1 output-table called: Subtotal

    tot_diff_Pol = 0.0            # total of the atomic polar contributions to solvation free enthalpy: wat - hex
    tot_diff_Apol = 0.0           # total of the atomic apolar contributions to solvation free enthalpy: wat - hex
    tot_diff_PolPlusApol = 0.0    # sum of polar and apolar: wat - hex



    for i in range(N):
           # alist[i][0] Atom number    (1 up to number of atoms in molecule)
           # alist[i][1] Chem. symbol
           # alist[i][2] CM2 chg.       (CM2 partial atomic charge)
           # alist[i][3] G_P (kcal/mol) (atomic polar contribution to solvation free energy)
           # alist[i][4] Area (Ang**2)  (surface area)
           # alist[i][5] Sigma (kcal/Ang**2)  (from AMSOL7.1 source code: Sigma = 1000D0*SRFACT(L)/ATAR(L))
           #                                  (i.e., Sigma[i] = 1000 * SS G_CDS[i]/Area[i], for atom i)
           #                                  ("Modeling free Energies of Solvation and Transfer", Computational Thermochemistry, Chapter 15,
           #                                  D.J. Giesen, D.G. Truhlar, 1998)
           # alist[i][6] SS G_CDS (kcal/mol) (atomic apolar contribution to solvation free energy)
           # alist[i][7] Subtotal (kcal/mol) (atomic solvation free energy (polar+apolar))
           # alist[i][8] M value


           # hexadecane:
           Chghex[i]  = float(atom_listhex[i][2])  # partial charges (CM2: Truhlar's charge model 2 charges) of each atom in hexadecane solvent
           Polhex[i]  = float(atom_listhex[i][3])  # atomic polar contributions to solvation free enthalpy in hexadecane: Polarization Free Energy (G_P) in kcal/mol
           SAA[i]     = float(atom_listhex[i][4])  # atomic surface area contibutions in Angstrom^2
           Sigmahex[i] = float(atom_listhex[i][5]) # Sigma (kcal/Ang**2)
           # in amsol-mod4 output::: 
           #Apolhex[i] = float(atom_listhex[i][5]); # atomic apolar contributions to solvation free enthalpy in hexadecane: SS G_CDS in kcal/mol 
           Apolhex[i] = float(atom_listhex[i][6])  # atomic apolar contributions to solvation free enthalpy in hexadecane: SS G_CDS in kcal/mol 

           # water:
           Polwat[i]  = float(atom_listwat[i][3])  # atomic polar contributions to solvation free enthalpy in water: Polarization Free Energy (G_P) in kcal/mol
           # in amsol-mod4 output::: 
           #Apolwat[i] = float(atom_listwat[i][5]) # atomic apolar contributions to solvation free enthalpy in water: SS G_CDS in kcal/mol
           Apolwat[i] = float(atom_listwat[i][6])  # atomic apolar contributions to solvation free enthalpy in water: SS G_CDS in kcal/mol

           ## sums over atomic contributions:

           Chghexsum += Chghex[i] # sum of the partial charges (CM2: Truhlar's charge model 2 charges) of each atom in hexadecane (!) solvent
           Polhexsum += Polhex[i] # sum of the atomic polar contributions to solvation free enthalpy: Polarization Free Energy (G_P)
           SAAsum += SAA[i]          # sum of the atomic surface area contibutions in Angstrom
                                     # is independent of the solvent under consideration
           # !!! former apolsum
           sum_Apolhex += Apolhex[i]    # sum of the atomic apolar contributions to solvation free enthalpy: SS G_CDS from hexadecane calculation (!!!!)
           # !!! former energy_sum
           # in amsol-mod4 output::: 
           #FreeEnergyHex_sum = FreeEnergyHex_sum + float(atom_listhex[i][6]); # sum of the atomic contributions to the total solvation free energy
                                                                               # for hexadecane (!!!) solvent: In AMSOL7.1 output-table called: Subtotal (in kcal/mol)
           FreeEnergyHex_sum += float(atom_listhex[i][7]); # sum of the atomic contributions to 
                                                           # for hexadecane (!!!) solvent: In AMSOL7.1 output-table called: Subtotal (in kcal/mol)
    #print tothex
    #print "len(tothex) :::", len(tothex)
    #print "tothex :::", tothex 
    #print "Chghexsum :::", Chghexsum
    #print "Polhexsum :::", Polhexsum
    #print "SAAsum :::", SAAsum
    #print "sum_Apolhex :::", sum_Apolhex
    #print "FreeEnergyHex_sum :::", FreeEnergyHex_sum
    #
    #
    # The information written in the "Total: ..." line close to the end of the hexadecane AMSOL7.1 output:
    #
    # tothex[0] = "Total:"
    # tothex[1] = CM2 chg.
    # tothex[2] = G_P, i.e. polar contribution to solvation free energy or also called Polarization free energy (in kcal/mol)
    # tothex[3] = Area (Ang**2) (in amsol-mod4: Area (CD) (Ang**2)): sum of all atomic contributions (!)
    # tothex[4] = SS G_CDS in kcal/mol (in amsol-mod4: G-CD)
    # tothex[5] = Subtotal in kcal/mol, i.e. sum of all atomic polar (G_P[i]) and apolar (SS G_CDS[i]) contributions to the solvation free energy
    #            Subtotal = solvation free energy
    #            Subtotal = "(5)  G-P-CDS(sol) = G-P(sol) + G-CDS(sol) = (2) + (4)             -XX.XXX kcal" 

    cs_coeff = (float(tothex[4]) - sum_Apolhex)/float(tothex[3])   # cs_coeff = total LS contribution / total AreaSAA
    #           (total SS G_CDS - sum of atomic all atomic apolar contributions to free energy of solvation)/ total Area
    #           = total CS/LS contribution divided by total Area

    #print ""
    #print "tothex[4] (SS G_CDS):::", float(tothex[4])
    #print "totwat[4] (SS G_CDS):::", float(totwat[4])
    #print "sum_Apolhex :::", sum_Apolhex
    #print "tothex[3] (Area) :::", float(tothex[3])
    #print "CS coeff ::: ", cs_coeff

    sum_csTimesSAA = 0.0

    # !!! WATER minus HEXADECANE: performing wat-hex difference !!! : 
    for i in range(N): # do substraction
           diff_Pol[i] = Polwat[i] - Polhex[i]                               # diff_Pol[i] is the difference of the atomic polar contribution to solvation free energy
                                                                             # obtained in water vs. the one obtained in hexadecane
                                                                             
           tot_diff_Pol = tot_diff_Pol + diff_Pol[i]                         # tot_diff_Pol is the sum of all differences of atomic polar contributions to solvation
                                                                             # free energy obtained in water vs. in hexadecane

           sum_csTimesSAA = sum_csTimesSAA + cs_coeff * SAA[i]               # sum_i cs_coeff * SAA[i] = sum_i (LS contribution/ total Area) * atomic Surface Area SAA[i] = LS contribution

                                                                             # - (Apolhex[i] + cs_coeff * SAA[i]): minus since: wat MINUS hex. cs_coeff and LS contribution only present in hex (Large Solvent)
           diff_Apol[i] = Apolwat[i] - (Apolhex[i] + cs_coeff * SAA[i])      # diff_Apol[i] is the difference of the atomic apolar contribution to solvation free energy
                                                                             # obtained in water vs. the one obtained in hexadecane minus cs_coeff * SAA[i]  
                                                                             
           tot_diff_Apol = tot_diff_Apol + diff_Apol[i]                      # tot_diff_Apol is sum of all atomic diff_Apol[i]
                                                                             
           diff_AtomicSolv[i] = diff_Pol[i] + diff_Apol[i]                   # diff_AtomicSolv[i] is the sum of the atomic diff_Pol[i] and diff_Apol[i] (see above for diff_Pol[i] and diff_Apol[i])
                                                                             # diff_AtomicSolv[i] considers the difference between water and hexadecane results for each atom!!!

           tot_diff_PolPlusApol = tot_diff_PolPlusApol + diff_AtomicSolv[i]  # tot_diff_PolPlusApol is the sum over all atomic contributions to the difference of solvation free energies
                                                                             # in water and hexadecane, i.e. the difference
                                                                             # of the final solvation free energies obtained in water and hexadecane 

    ## write out the solvation file 
    file = open(outputprefix+'.solv','w')
    file.write( "%s %3d %4.1f %8.2f %8.2f %8.2f %8.2f\n" % (name,numatom,float(totwat[1]),tot_diff_Pol,float(totwat[3]),tot_diff_Apol,tot_diff_PolPlusApol))  # formal charge is the same for water or hexadecane
    for i in range(N):
           file.write("%8.4f%8.2f%7.2f%8.2f%8.2f\n" % (Chghex[i],diff_Pol[i],SAA[i],diff_Apol[i],diff_AtomicSolv[i]))
    file.close()

    if VERBOSE:
        print ("")
        print ("**** The function diff_amsol71_files() was finished. ****")
        print ("")

    return

#################################################################################################################
#################################################################################################################

def modify_charges_mol2_file(mol2file, atom_list_hex, outputprefix):
    ## read in mol2 file

    if VERBOSE:
        print ("")
        print ("**** starting the function modify_charges_mol2_file() *****")
        print ("")
        print ("     CM2 charges from an AMSOL7.1 SM5.42R calculation in hexadecane (!!!) ")
        print ("     are written to mol2-file. Former charges are overwritten.")

    #print ("I AM HERE %s"%mol2file)
    mol = mol2amsol.read_Mol2_file(mol2file)[0] 

    n = len(atom_list_hex)
    if n != len(mol.atom_list):
       print ("Error: n != len(mol.atom_list) : " + str(n) + " !=" + str(len(mol.atom_list)))
       exit()

    
    for i in range(n):
        charge = atom_list_hex[i][2]
        if VERBOSE:
            print (mol.atom_list[i].Q, charge)
        mol.atom_list[i].Q = float(charge)

    filename = outputprefix + '.mol2'
    mol2amsol.write_mol2(mol,filename)
    
    if VERBOSE:
        print ("")
        print ("**** The function modify_charges_mol2_file() was finished. ****")
        print ("")

    return
    
#################################################################################################################
#################################################################################################################
def main():

    if VERBOSE:
        print ("")
        print ("**** entering the main program in process_amsol71_mol2.py ****")
        print ("")

    if len(sys.argv) != 5: # if no input
        print (" This script needs the following:")
        print (" (1) amsol water output file ")
        print (" (2) amsol hex output file ")
        print (" (3) amsol hex output file ")
        print (" (4) outputprefix ")
        return

    filenamewat    = sys.argv[1]
    filenamehex    = sys.argv[2]
    mol2file       = sys.argv[3]
    outputprefix   = sys.argv[4]

    atom_list_wat,tot_wat,name_wat,numat_wat = process_amsol_file(filenamewat,outputprefix,"wat")
    atom_list_hex,tot_hex,name_hex,numat_hex = process_amsol_file(filenamehex,outputprefix,"hex")
    

    # tot_wat and tot_hex are lists:
    # ( see def process_amsol_file(...) above)
    #
    # tot_wat[0] or tot_hex[0] : total_line[0] = "Total:"
    # tot_wat[1] or tot_hex[1] : total_line[1] = CM2 chg.
    # tot_wat[2] or tot_hex[2] : total_line[2] = G_P, i.e. polar contribution to solvation free energy or also called Polarization free energy (in kcal/mol)
    # tot_wat[3] or tot_hex[3] : total_line[3] = Area (Ang**2) (in amsol-mod4: Area (CD) (Ang**2))
    # tot_wat[4] or tot_hex[4] : total_line[4] = SS G_CDS in kcal/mol (in amsol-mod4: G-CD)
    # tot_wat[5] or tot_hex[5] : total_line[5] = Subtotal in kcal/mol, i.e. sum of all atomic polar (G_P[i]) and apolar (SS G_CDS[i]) contributions to the solvation free energy
    #                            Subtotal = solvation free energy
    #                            Subtotal = "(5)  G-P-CDS(sol) = G-P(sol) + G-CDS(sol) = (2) + (4)             -XX.XXX kcal"
    


    if (name_hex != name_wat or numat_hex != numat_wat):
        print ("Error: Name or Atom counts do not agree")
    elif VERBOSE:
        print ("Wat-name      = " + name1 )
        print ("Wat-atom cout = " + str(numat1))
        print ("Hex-name      = " + name2 )
        print ("Hex-atom cout = " + str(numat2))

    if VERBOSE:
        print ("")
        print ("just before diff_amsol71_files() function")
        print ("")
    diff_amsol71_files(atom_list_wat, tot_wat, atom_list_hex, tot_hex, name_wat, numat_wat, outputprefix)

    if VERBOSE:
        print ("")
        print ("just before modify_charges_mol2_file() function")
        print ("")
    modify_charges_mol2_file(mol2file, atom_list_hex, outputprefix) 

    if VERBOSE:
        print ("")
        print ("**** The main program in process_amsol71_mol2.py was finished for ****")
        print ("%s and %s." % (filenamewat, filenamehex))
        print ("*******************************************************************")
        print ("")
    #
    return 
#################################################################################################################
#################################################################################################################
main()
