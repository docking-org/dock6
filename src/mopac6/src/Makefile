# Makefile for mopac

include ../../../install/rules.h
include ../../../install/config.h

PROGRAM = mopac$(SFX)
BIN     = ../../../bin

OBJS    = aababc.o  readmo.o writmo.o mopac.o analyt.o \
anavib.o  axis.o  block.o  bonds.o  brlzon.o  calpar.o \
capcor.o  cdiag.o  chrge.o  cnvg.o  compfg.o  datin.o \
dcart.o  delmol.o  delri.o  denrot.o  densit.o  depvar.o \
deri0.o  deri1.o  deri2.o  deri21.o  deri22.o  deri23.o \
deritr.o  deriv.o  dernvo.o  ders.o  dfock2.o  dfpsav.o \
diag.o  diat.o  diat2.o  diis.o  dijkl1.o  dijkl2.o \
dipind.o  dipole.o  dofs.o  dot.o  drc.o  drcout.o \
ef.o  enpart.o  exchng.o  ffhpol.o  flepo.o  fmat.o \
fock1.o  fock2.o  force.o  formxy.o  forsav.o  frame.o \
freqcy.o  geout.o  geoutg.o  getgeg.o  getgeo.o  getsym.o \
gettxt.o  gmetry.o  gover.o  grid.o  h1elec.o  haddon.o \
hcore.o  helect.o  hqrii.o  ijkl.o  interp.o  iter.o \
jcarin.o  linmin.o  local.o  locmin.o  mamult.o  matout.o \
matpak.o  meci.o  mecid.o  mecih.o  mecip.o  moldat.o \
molval.o  mullik.o  mult.o  nllsq.o  nuchar.o  parsav.o \
partxy.o  pathk.o  paths.o  perm.o  polar.o  powsav.o \
powsq.o  prtdrc.o  quadr.o  react1.o  reada.o  refer.o \
repp.o  rotat.o  rotate.o  rsp.o  search.o  second.o \
setupg.o  solrot.o  swap.o  symtry.o  thermo.o  timer.o \
update.o  vecprt.o  wrtkey.o  wrttxt.o  xyzint.o 

SRCS     = aababc.f  analyt.f  anavib.f  axis.f  block.f \
bonds.f  calpar.f  capcor.f  chrge.f  cnvg.f  compfg.f \
dcart.f  delmol.f  delri.f  denrot.f  depvar.f  deriv.f \
ders.f  dfpsav.f  diat.f  diat2.f  dipind.f  dipole.f \
dot.f  drc.f  drcout.f  enpart.f  exchng.f  ffhpol.f \
flepo.f  fmat.f  fock1.f  force.f  forsav.f  frame.f \
freqcy.f  geout.f  getgeo.f  getsym.f  gmetry.f  gover.f \
grid.f  h1elec.f  haddon.f  hcore.f  helect.f  hqrii.f \
ijkl.f  interp.f  iter.f  linmin.f  local.f  locmin.f \
mamult.f  matout.f  meci.f  moldat.f  molval.f  mullik.f \
mult.f  nllsq.f  nuchar.f  parsav.f  paths.f  perm.f \
polar.f  powsav.f  powsq.f  prtdrc.f  quadr.f  react1.f \
reada.f  refer.f  repp.f  rotat.f  search.f  second.f \
setupg.f  solrot.f  swap.f  symtry.f  thermo.f  update.f \
vecprt.f  wrtkey.f  xyzint.f

all:		$(PROGRAM)

$(PROGRAM):     $(OBJS)
		$(FC) $(LDFLAGS) $(OBJS) $(LIBS) -o $(PROGRAM)

install:	$(PROGRAM)
		mv $(PROGRAM) $(BIN)

clean:
	/bin/rm -f $(OBJS) 

realclean: clean
	/bin/rm -f $(PROGRAM) 

uninstall:
	cd $(BIN) && /bin/rm -f $(PROGRAM)

