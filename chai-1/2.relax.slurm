#!/bin/bash
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

# D0C is optional. Resolve its project-local chemistry files before changing
# into the sample directory, but load them only for inputs that contain D0C.
d0c_param_flags=()
if awk '
    ($1 == "ATOM" || $1 == "HETATM") && $6 == "D0C" { found=1; exit }
    END { exit !found }
' "$sample_dir"/pred.model_idx_*.cif; then
    d0c_params=$PWD/rosetta_params/D0C.params
    d0c_bonded_params=$PWD/rosetta_params/bondlength_bondangle_d0c
    rosetta_database=$ROSETTAROOT/main/database

    for required_file in \
        "$d0c_params" \
        "$d0c_bonded_params/default-improper.txt" \
        "$rosetta_database/chemical/residue_type_sets/fa_standard/residue_types.txt"; do
        if [ ! -r "$required_file" ]; then
            echo "Required D0C Rosetta input is missing or unreadable: $required_file" >&2
            exit 1
        fi
    done

    d0c_param_flags=(
        -extra_res_fa "$d0c_params"
        -in:path:database "$PWD/rosetta_params" "$rosetta_database"
        -score:bonded_params_dir bondlength_bondangle_d0c
    )
fi

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
 "${d0c_param_flags[@]}" \
 >/dev/null \
 2>err.log
# -out:path:all relax \

fi
