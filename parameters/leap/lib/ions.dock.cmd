# DOCK Amber Score Ion Library
#  Scott Brozell August 2008

# Create residues and units for some ions not treated by the Amber 
#  force fields.  This is based on ions94.
#  Some of these ions are already treated by DOCK:
#  these are alkaline earth divalent, monoatomic ions for Sr;
#  sensibly, the names chosen, SR, correspond to the
#  common names of these ions in pdb's.
#  Some of these ions are not yet treated by DOCK:
#  Hg(II), Mn(II)

clearVariables
logFile ions.dock.log

i = createAtom   SR  Sr  2.0
set i    element Sr
set i    position { 0 0 0 }
r = createResidue SR
add r i
SR = createUnit SR
add SR r
saveOff SR ./ions.dock.lib


i = createAtom   HG  Hg  2.0
set i    element Hg
set i    position { 0 0 0 }
r = createResidue HG
add r i
HG = createUnit HG
add HG r
saveOff HG ./ions.dock.lib

i = createAtom   MN  Mn  2.0
set i    element Mn
set i    position { 0 0 0 }
r = createResidue MN
add r i
MN = createUnit MN
add MN r
saveOff MN ./ions.dock.lib

quit
