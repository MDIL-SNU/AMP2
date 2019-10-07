echo 'node information' > nodefile
NPROC=16

### set path of config.yaml ###
conf=./config.yaml
###############################

### Do not change #############
src_path=`grep 'src_path' $conf | tr -s ' ' | cut -d " " -f 3`
###############################

python $src_path/main.py $conf nodefile $NPROC >& stdout.x
