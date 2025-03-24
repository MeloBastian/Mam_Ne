#!/usr/bin/env bash
# run_busco.sbatch --------------------------------------------------------

#SBATCH --constraint="skylake|haswell|broadwell"
#SBATCH --job-name="busco_${1}${2}"
#SBATCH --exclude=pbil-deb[6-20,38],lbbe-t4cloud
#SBATCH -N1
#SBATCH -n1
#SBATCH -c8
#SBATCH --mem 32G
#SBATCH --time "60:00:00"

set -x

outdir="busco_${1}_${2}"
localdir="busco_analysis${1}_${2}/"
cwd=${PWD}
metaeuk="${outdir}/run_mammalia${2}_odb10/metaeuk_output"

cleanup() { 
 # if [ ! -f "${outdir}/DONE.flag" ]; then #always verified
    rm -rf "${metaeuk}/tmp" 
    rm -rf "${metaeuk}/refseq_db_rerun.faa"
    echo "metaeuk temp file removed"
    rsync -avhr ${outdir} ${cwd}/busco_/
    cd ..
    rm -rf ${localdir}
 # fi
}

trap cleanup EXIT

export PATH="/beegfs/data/soft/singularity-3.1.1/bin:$PATH"
cpus=8
assembly="${1}"

cd "/data/mbastian/"
mkdir ${localdir}
mkdir ${localdir}/busco_downloads
rsync -avhr ${cwd}/busco_downloads/lineages/mammalia${2}_odb10 ${localdir}/busco_downloads/lineages/
rsync -avhr ${cwd}/${assembly} ${localdir}
cd ${localdir} 

singularity exec https://depot.galaxyproject.org/singularity/busco:5.0.0--py_1 busco\
  --offline \
  --in ${assembly} \
  --out ${outdir} \
  --mode genome \
  --lineage mammalia${2}_odb10 \
  --cpu $cpus \
  --force \
  

