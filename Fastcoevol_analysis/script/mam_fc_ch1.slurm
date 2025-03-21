#!/bin/bash
#SBATCH --time=14:00:00
#SBATCH --nodes=1
#SBATCH --job-name="mamfc1"
#SBATCH --ntasks-per-node=1
#SBATCH --mem=1gb
#SBATCH -o mam_fc_ch1.out
#SBATCH -e mam_fc_ch1.err

../../../../bayescode/data/fastcoevol -ss mam_empiricaldnds_merge_dsomsuffstat.dsompathsuffstat -c thvwithoutgenerationtime -r rootfile -t 1007forphylo.conc.treefile_rooted -wn mam_fc_withoutgenerationtime_6002maskedgenes_ch1
