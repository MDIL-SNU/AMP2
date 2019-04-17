#!/bin/sh
#PBS -l nodes=x01:ppn=16
#PBS -N AMP2_test
##PBS -l walltime=48:00:00

cd $PBS_O_WORKDIR
cat $PBS_NODEFILE > nodefile
NPROC=`wc -l < $PBS_NODEFILE`

### set path of config.yaml ###
conf=./config.yaml
###############################

### Do not change #############
src_path=`grep 'src_path' $conf | tr -s ' ' | cut -d " " -f 3`
###############################

python $src_path/main.py $conf nodefile $NPROC >& stdout.x
