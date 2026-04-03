#!/bin/sh 
#SBATCH -p cpu
#SBATCH -J relax
#SBATCH -o relax.log
#SBATCH -e relax.err
#SBATCH -a 1
#SBATCH --mem 2g
##SBATCH -x node12

# set -e
. ~/.bashrc

sample_dir=`sed -n ${SLURM_ARRAY_TASK_ID}p outputs.list`
dir=$(dirname $(dirname $sample_dir))
xml=$PWD/relax.xml

# check if D-chains come at the front or the end
n_seqs=`cat $dir/example.fasta| grep '>'| wc -l`
n_d_seqs=`cat $dir/example.fasta| grep ')(D'| wc -l`
d_seqs_front=`head -n 2 $dir/example.fasta|grep ')(D' |wc -l`
if [ $d_seqs_front -gt 0 ];then
    ((start=n_d_seqs+1))  # front
else
    ((start=n_seqs-n_d_seqs+1))
fi

chain_grp2_ids=""
for i in `seq $start $n_seqs`;do
    chain_grp2_ids+=",$i"
done

cd $sample_dir
prefix=`basename $sample_dir`
pdbs=(`ls ${prefix}_pred.model_idx_*_0001.pdb`)
if [ ${#pdbs[*]} -ne 5 ];then

rm -f ${prefix}_score.sc
ls *.cif >relax.list

$rosettabin/rosetta_scripts$rosettasuffix \
 -parser:protocol $xml \
 -l relax.list \
 -script_vars chain_grp2_ids=${chain_grp2_ids:1} \
 -nstruct 1 \
 -ex1 \
 -ex2aro \
 -out:prefix ${prefix}_ \
 -overwrite \
 -beta \
 >/dev/null \
 2>err.log
# -out:path:all relax \

fi
