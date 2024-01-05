#!/bin/bash

gcc LuisArmando-Projeto-Sequencial.c -o LuisArmando-Projeto-Sequencial -lm
gcc -O2 -ffinite-math-only -mtune=native -funroll-loops -foptimize-register-move -m64 -fopenmp LuisArmando-Projeto-Otimizado.c -o LuisArmando-Projeto-Otimizado -lm
gcc -O2 -ffinite-math-only -mtune=native -funroll-loops -foptimize-register-move -m64 -fopenmp LuisArmando-Projeto-OpenMP.c -o LuisArmando-Projeto-OpenMP -lm

echo "Compilação Finalizada"


