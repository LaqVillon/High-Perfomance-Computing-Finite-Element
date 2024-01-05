#!/bin/bash

export FI_SOCKETS_IFACE=eth0
export FI_PROVIDER=tcp
ulimit -s unlimited
ulimit -v unlimited
# mpicc -O2 -ffinite-math-only -mtune=native -funroll-loops -foptimize-register-move -m64 LuisArmando-Programa-MPI.c -o LuisArmando-Programa-MPI -lm #
# mpicc -O2 -ffinite-math-only -mtune=native -funroll-loops -foptimize-register-move -m64 LuisArmando-Programa-Sequencial.c -o LuisArmando-Programa-Sequencial -lm #
mpiicc -ipo -xhost -ffinite-math-only LuisArmando-Programa-MPI.c -o LuisArmando-Programa-MPI -lm 
mpiicc -ipo -xhost -ffinite-math-only LuisArmando-Programa-Sequencial.c -o LuisArmando-Programa-Sequencial -lm

echo "Compilação Finalizada"


