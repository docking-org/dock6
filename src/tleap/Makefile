#  This is a makefile to just make text-related programs, primarily teLeap
#
include ../../install/config.h

BIN = ../../bin

.SUFFIXES:
.SUFFIXES:  .c .o

.c.o: $<
	$(CC) -c $(CFLAGS) -DSYSV $<

#
#  Absolutely necessary routines
#

OA1= basics.o sysdepend.o stringExtra.o varArray.o  getline.o avl.o
OA = $(OA1)

#
#  General tools
#

OG1= pdb_format.o pdb_read.o pdb_sprntf.o
OG2= pdb_sscanf.o pdb_write.o
OG3= vector.o zMatrix.o sort.o bag.o hash.o
OG4= dictionary.o database.o nVector.o ring.o matrix.o
OG5= fortran.o displayer.o
OG = $(OG1) $(OG2) $(OG3) $(OG4) $(OG5)

#
#  Object hierarchy objekts
#

OO1= assoc.o atom.o byteArray.o
OO2= collection.o container.o internal.o
OO3= list.o loop.o molecule.o
OO4= oDouble.o oInteger.o oString.o
OO5= objekt.o parmSet.o residue.o unit.o unitio.o
OO6= graphUtil.o select.o
OO = $(OO1) $(OO2) $(OO3) $(OO4) $(OO5) $(OO6)

#
#  leap-specific objects and tools
#

OL1= amber.o build.o elements.o library.o chirality.o
OL2= minimizer.o model.o parmLib.o pdbFile.o tools.o
OL3= variables.o parser.o help.o helptext.o octree.o
OL4= commands.o mathop.o block.o restraint.o hybrid.o tripos.o
OL = $(OL1) $(OL2) $(OL3) $(OL4)

#
#  Terminal interface specifics
#

OC1= tLeap.o
OC  = $(OC1)

TLEAP_OBJ= $(OA) $(OG) $(OO) $(OL) $(OC)
teLeap$(SFX): $(TLEAP_OBJ)
	$(CC) -o $@ $(TLEAP_OBJ) -lm

#  material to make a text-only tar file:

SRC= basics.c sysdepend.c stringExtra.c varArray.c  getline.c avl.c \
	pdb_format.c pdb_read.c pdb_sprntf.c \
	pdb_sscanf.c pdb_write.c \
	vector.c zMatrix.c sort.c bag.c hash.c \
	dictionary.c database.c nVector.c ring.c matrix.c \
	fortran.c displayer.c \
	assoc.c atom.c byteArray.c \
	collection.c container.c internal.c \
	list.c loop.c molecule.c \
	oDouble.c oInteger.c oString.c \
	objekt.c parmSet.c residue.c unit.c unitio.c \
	graphUtil.c select.c tripos.c \
	amber.c build.c elements.c library.c chirality.c \
	minimizer.c model.c parmLib.c pdbFile.c tools.c \
	variables.c help.c helptext.text octree.c \
	commands.c mathop.c block.c restraint.c hybrid.c \
	tLeap.c parser.y \
	amber.h assoc.h atom.h avl.h avl2.h bag.h \
	basics.h block.h build.h byteArray.h chirality.h classes.h \
	collection.h commands.h container.h database.h dictionary.h displayer.h \
	elements.h fortran.h function.h getline.h graphUtil.h hash.h \
	help.h hybrid.h internal.h leap.h library.h list.h \
	loop.h machine.h mathop.h matrix.h mess.h minimizer.h \
	model.h molecule.h nVector.h oDouble.h oInteger.h oString.h \
	objekt.h octree.h parmLib.h parmSet.h parser.h pdb.h \
	pdbFile.h pdb_int.h residue.h restraint.h rgb2hls.h ring.h \
	select.h sort.h stringExtra.h threed.h tools.h tripos.h \
	unit.h unitio.h varArray.h vector.h \
	Makefile.tleap

tarfile: 
	tar cvf tleap.tar $(SRC)

# end of tarfile section

OUMH= utilMakeHelp.o
UMH_OBJ = $(OA) $(OG) $(OUMH)
utilMakeHelp: $(UMH_OBJ)
	$(CC) -o $@ $(UMH_OBJ) -lm

parser.c:       parser.y parser.h
	$(YACC) parser.y
	@echo "(The warning 'conflicts: 6 shift/reduce' is normal)"
	-rm -f parser.c
	mv y.tab.c parser.c

helptext.c: helptext.text utilMakeHelp
	./utilMakeHelp helptext.text helptext.c

amber.o: basics.h
amber.o: vector.h
amber.o: classes.h
amber.o: leap.h
amber.o: tools.h
amber.o: amber.h
assoc.o: basics.h
assoc.o: classes.h
atom.o: basics.h
atom.o: classes.h
atom.o: elements.h
avl.o: avl.h
bag.o: basics.h
bag.o: bag.h
basics.o: basics.h
basics.o: varArray.h
block.o: basics.h
block.o: block.h
build.o: basics.h
build.o: vector.h
build.o: matrix.h
build.o: zMatrix.h
build.o: classes.h
build.o: chirality.h
build.o: minimizer.h
build.o: parmLib.h
build.o: model.h
build.o: build.h
byteArray.o: basics.h
byteArray.o: classes.h
chirality.o: basics.h
chirality.o: classes.h
chirality.o: chirality.h
collection.o: basics.h
collection.o: classes.h
commands.o: basics.h
commands.o: vector.h
commands.o: matrix.h
commands.o: classes.h
commands.o: dictionary.h
commands.o: database.h
commands.o: library.h
commands.o: parmLib.h
commands.o: pdbFile.h
commands.o: help.h
commands.o: parser.h
commands.o: tools.h
commands.o: amber.h
commands.o: commands.h
commands.o: leap.h
commands.o: octree.h
container.o: basics.h
container.o: vector.h
container.o: classes.h
container.o: container.h
database.o: basics.h
database.o: dictionary.h
database.o: database.h
dictionary.o: basics.h
dictionary.o: hash.h
dictionary.o: dictionary.h
displayer.o: basics.h
displayer.o: displayer.h
elements.o: basics.h
elements.o: elements.h
fortran.o: basics.h
function.o: basics.h
function.o: varArray.h
function.o: function.h
getline.o: getline.h
graphUtil.o: basics.h
graphUtil.o: classes.h
graphUtil.o: graphUtil.h
hash.o: basics.h
hash.o: hash.h
help.o: basics.h
help.o: dictionary.h
help.o: help.h
hybrid.o: basics.h
internal.o: basics.h
internal.o: classes.h
internal.o: ring.h
library.o: basics.h
library.o: classes.h
library.o: dictionary.h
library.o: database.h
library.o: parmLib.h
library.o: library.h
list.o: basics.h
list.o: classes.h
loop.o: basics.h
loop.o: classes.h
mathop.o: basics.h
mathop.o: vector.h
mathop.o: matrix.h
mathop.o: classes.h
matrix.o: basics.h
matrix.o: vector.h
matrix.o: matrix.h
mess.o: basics.h
mess.o: function.h
mess.o: mess.h
minimizer.o: basics.h
minimizer.o: varArray.h
minimizer.o: nVector.h
minimizer.o: classes.h
minimizer.o: minimizer.h
model.o: basics.h
model.o: classes.h
model.o: chirality.h
model.o: model.h
molecule.o: basics.h
molecule.o: classes.h
nVector.o: basics.h
nVector.o: nVector.h
oDouble.o: basics.h
oDouble.o: classes.h
oInteger.o: basics.h
oInteger.o: classes.h
oString.o: basics.h
oString.o: classes.h
objekt.o: basics.h
objekt.o: classes.h
octree.o: classes.h
octree.o: octree.h
parmLib.o: basics.h
parmLib.o: classes.h
parmLib.o: parmLib.h
parmSet.o: basics.h
parmSet.o: classes.h
parmSet.o: dictionary.h
parmSet.o: database.h
parmSet.o: amber.h
pdbFile.o: basics.h
pdbFile.o: classes.h
pdbFile.o: dictionary.h
pdbFile.o: database.h
pdbFile.o: pdb.h
pdbFile.o: varArray.h
pdbFile.o: avl.h
pdbFile.o: sort.h
pdbFile.o: matrix.h
pdbFile.o: model.h
pdbFile.o: minimizer.h
pdbFile.o: leap.h
pdb_format.o: pdb_int.h
pdb_read.o: basics.h
pdb_read.o: pdb_int.h
pdb_write.o: pdb_int.h
residue.o: basics.h
residue.o: classes.h
restraint.o: basics.h
restraint.o: classes.h
restraint.o: restraint.h
ring.o: basics.h
ring.o: ring.h
select.o: basics.h
select.o: classes.h
select.o: minimizer.h
select.o: select.h
sort.o: basics.h
stringExtra.o: basics.h
sysdepend.o: basics.h
tLeap.o: basics.h
tLeap.o: classes.h
tLeap.o: dictionary.h
tLeap.o: parser.h
tLeap.o: block.h
tLeap.o: leap.h
tLeap.o: getline.h
test_dataBase.o: basics.h
test_dataBase.o: database.h
test_parmSet.o: basics.h
test_parmSet.o: classes.h
test_parmSet.o: database.h
test_parmSet.o: amber.h
test_varArray.o: basics.h
test_varArray.o: varArray.h
tools.o: basics.h
tools.o: vector.h
tools.o: matrix.h
tools.o: classes.h
tools.o: dictionary.h
tools.o: database.h
tools.o: parser.h
tools.o: leap.h
tools.o: tools.h
tripos.o: basics.h
tripos.o: classes.h
tripos.o: tripos.h
unit.o: basics.h
unit.o: vector.h
unit.o: classes.h
unit.o: restraint.h
unit.o: bag.h
unit.o: dictionary.h
unit.o: database.h
unit.o: parmLib.h
unit.o: unitio.h
unitio.o: basics.h
unitio.o: vector.h
unitio.o: classes.h
unitio.o: restraint.h
unitio.o: bag.h
unitio.o: dictionary.h
unitio.o: database.h
unitio.o: parmLib.h
unitio.o: avl.h
varArray.o: basics.h
varArray.o: varArray.h
variables.o: basics.h
variables.o: classes.h
variables.o: dictionary.h
variables.o: leap.h
vector.o: basics.h
vector.o: vector.h
zMatrix.o: basics.h
zMatrix.o: vector.h
zMatrix.o: matrix.h
zMatrix.o: zMatrix.h

install:  teLeap$(SFX)
	mv teLeap$(SFX) $(BIN)
	cd ../..; DOCKHOME=`pwd`;
	sed "s?DOCKHOME?$(DOCKHOME)?" tleap > $(BIN)/tleap
	chmod +x $(BIN)/tleap

clean::
	/bin/rm -f teLeap$(SFX) utilMakeHelp$(SFX) *.o parser.c helptext.c

uninstall::
	cd $(BIN) && /bin/rm -f teLeap$(SFX) tleap

