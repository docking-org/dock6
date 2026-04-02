#!/bin/bash

case "$1" in
  "am1bcc" | "amberize_complex" | "amberize_ligand" | "amberize_receptor" | "antechamber" | "atomtype" | "bondtype" | \
  "chemgrid" | "distmap" | "espgen" | "grid" | "grid-convrds" | "grid-convert" | \
  "make_phimap" | "mopac" | "mopac.sh" | \
  "nchemgrid_GB" | "nchemgrid_SA" | "parmcal" | "parmchk" | "prepgen" | \
  "resp" | "respgen" | "sevsolv" | "showbox" | "showsphere" | "sphere_selector" | "sphgen" | "solvgrid" | "solvmap" | \
  "teLeap" | "tleap" | "dock6" | "dock6.rdkit" | "dock6.mpi" | "mpirun")
    exec "$1" "${@:2}"
    ;;
  *)
    echo "$1 is an unknown binary..." && \
    echo "###################" && \
    echo "-----Available Binaries-----:" && \
    echo "$(ls -I *.docker -I *py -I *pyc -I *pl -I *sh /app/dock6/bin)"  && \
    echo "mpirun"  && \
    echo "-----Example mpi runs-----" && \
    echo "\"mpirun -n 10 dock6.mpi -i dock.in -o dock.out\""
    exit 1
    ;;
esac
