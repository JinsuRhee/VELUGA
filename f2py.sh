##-----
## FILE LIST FOR FORTRAN
##-----
FILES='get_ptcl_py.f90 get_cell_py.f90 get_amr_py.f90 find_domain_py.f90'

##-----
## WRAPPER
##-----
F2PY=f2py

#export CFLAGS="-fpermissive"
##-----
## FORTRAN COMPILER
##
## 	-- TODO --
##	How to specify this?
##-----
FORT=/home/jinsu/anaconda3/envs/rur/bin/x86_64-conda-linux-gnu-gfortran

##-----
## MACHINE SETTINGS?
##-----
MACHINE="gc"

##-----
## GO
##-----
BASEDIR=$(dirname "$0")
cd $BASEDIR/src/fortran

##----- Compile interanl modules first
$FORT -fopenmp -mcmodel=large -O2 -c -fPIC read_ramses_py.f90 -o read_ramses_py.o

##----- COMPILE WITH WRAPPER
for f in $FILES
do
  bn=$(basename "$f" .f90)
  echo -e "\n\n\nCompiling $bn\n\n\n"
  $F2PY -m $bn --f90flags='-fopenmp -O2' -I read_ramses_py.o -lgomp -c $f  
done
rm *.o

##-----
## crash happens when env gfortran != conda gfortran
##
## default - using gfortran library on your conda env


