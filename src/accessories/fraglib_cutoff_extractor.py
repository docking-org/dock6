#Editor: John Bickel
#Purpose: Keeps only fragments from the fragment library based on a given cutoff (13,000 default)
#Edit date: 03/04/2021

#Returns a list of all the lines of the molecules in a given multi-mol2
def extract_all_molecules(filename):
    line_list=[]
    molecule_list=[]
    with open(filename,"r") as f:
        append_val=0
        for line in f:
            #TYPE: is the first line in each header
            if "TYPE:" in line:
                append_val+=1 #shoddy way of tracking if we're at the second molecule
                if append_val>=2: #if we are, start appending things
                    append_next=True
                    molecule_list.append(line_list)
                    
                line_list=[]
                line_list.append(line)
            else:
                line_list.append(line)
    return(molecule_list) #returns a list

#takes all the fragments and only keeps the ones > freq_cutoff
def remake_by_cutoff(filename,outfile,freq_cutoff):
    #pulls all fragments in the file into a list
    fragments=extract_all_molecules(filename)
    with open(outfile,"w") as of:
        for frag in fragments:
            #frag[2] is the FREQ: header line with basic output
            if int(frag[2].split(":")[1].strip()) > freq_cutoff:
                for line in frag:
                    #writes to line if above freq_cutoff
                    of.write(line)
    

#Raw fragment library outputs
lnk_file="fraglib_linker.mol2"
scf_file="fraglib_scaffold.mol2"
sid_file="fraglib_sidechain.mol2"

#Outfile names of your choice
lnk_outfile="fraglib_linker_cutoff.mol2"
scf_outfile="fraglib_scaffold_cutoff.mol2"
sid_outfile="fraglib_sidechain_cutoff.mol2"

#Whatever you want the cutoff to be
freq_cutoff=13000

remake_by_cutoff(lnk_file,lnk_outfile,freq_cutoff)
remake_by_cutoff(scf_file,scf_outfile,freq_cutoff)
remake_by_cutoff(sid_file,sid_outfile,freq_cutoff)
