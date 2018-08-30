#!/bin/bash
#PBS -N salmon_seed -l walltime=48:00:00 -t 0-35
DIR='/work/maa146/EffectiveLength/salmon-edgeR/'
unset libs
declare -A libs
libs['A2-10-R1']='raw/seed_for_eflen_paper/A2-10-R1.cut.fq.gz'
libs['A2-10-R2']='raw/seed_for_eflen_paper/A2-10-R2.cut.fq.gz'
libs['A2-10-R3']='raw/seed_for_eflen_paper/A2-10-R3.cut.fq.gz'
libs['A2-20-R1']='raw/seed_for_eflen_paper/A2-20-R1.cut.fq.gz'
libs['A2-20-R2']='raw/seed_for_eflen_paper/A2-20-R2.cut.fq.gz'
libs['A2-20-R3']='raw/seed_for_eflen_paper/A2-20-R3.cut.fq.gz'
libs['A2-30-R1']='raw/seed_for_eflen_paper/A2-30-R1.cut.fq.gz'
libs['A2-30-R2']='raw/seed_for_eflen_paper/A2-30-R2.cut.fq.gz'
libs['A2-30-R3']='raw/seed_for_eflen_paper/A2-30-R3.cut.fq.gz'
libs['A2-40-R1']='raw/seed_for_eflen_paper/A2-40-R1.cut.fq.gz'
libs['A2-40-R2']='raw/seed_for_eflen_paper/A2-40-R2.cut.fq.gz'
libs['A2-40-R3']='raw/seed_for_eflen_paper/A2-40-R3.cut.fq.gz'
libs['ADs-10-R1']='raw/seed_for_eflen_paper/ADs-10-R1.cut.fq.gz'
libs['ADs-10-R2']='raw/seed_for_eflen_paper/ADs-10-R2.cut.fq.gz'
libs['ADs-10-R3']='raw/seed_for_eflen_paper/ADs-10-R3.cut.fq.gz'
libs['ADs-20-R1']='raw/seed_for_eflen_paper/ADs-20-R1.cut.fq.gz'
libs['ADs-20-R3']='raw/seed_for_eflen_paper/ADs-20-R3.cut.fq.gz'
libs['ADs-30-R1']='raw/seed_for_eflen_paper/ADs-30-R1.cut.fq.gz'
libs['ADs-30-R2']='raw/seed_for_eflen_paper/ADs-30-R2.cut.fq.gz'
libs['ADs-30-R3']='raw/seed_for_eflen_paper/ADs-30-R3.cut.fq.gz'
libs['ADs-40-R1']='raw/seed_for_eflen_paper/ADs-40-R1.cut.fq.gz'
libs['ADs-40-R2']='raw/seed_for_eflen_paper/ADs-40-R2.cut.fq.gz'
libs['ADs-40-R3']='raw/seed_for_eflen_paper/ADs-40-R3.cut.fq.gz'
libs['D5-10-R1']='raw/seed_for_eflen_paper/D5-10-R1.cut.fq.gz'
libs['D5-10-R2']='raw/seed_for_eflen_paper/D5-10-R2.cut.fq.gz'
libs['D5-10-R3']='raw/seed_for_eflen_paper/D5-10-R3.cut.fq.gz'
libs['D5-20-R1']='raw/seed_for_eflen_paper/D5-20-R1.cut.fq.gz'
libs['D5-20-R2']='raw/seed_for_eflen_paper/D5-20-R2.cut.fq.gz'
libs['D5-20-R3']='raw/seed_for_eflen_paper/D5-20-R3.cut.fq.gz'
libs['D5-30-R1']='raw/seed_for_eflen_paper/D5-30-R1.cut.fq.gz'
libs['D5-30-R2']='raw/seed_for_eflen_paper/D5-30-R2.cut.fq.gz'
libs['D5-30-R3']='raw/seed_for_eflen_paper/D5-30-R3.cut.fq.gz'
libs['D5-40-R1']='raw/seed_for_eflen_paper/D5-40-R1.cut.fq.gz'
libs['D5-40-R2']='raw/seed_for_eflen_paper/D5-40-R2.cut.fq.gz'
libs['D5-40-R3']='raw/seed_for_eflen_paper/D5-40-R3.cut.fq.gz'
cd $DIR/seed
ROOT=$(git rev-parse --show-toplevel)

PATH=$PATH:$DIR/Salmon-latest_linux_x86_64/bin

KEYS=("${!libs[@]}")
name="${KEYS[$PBS_ARRAYID]}"

salmon quant --no-version-check \
             -p $PBS_NUM_PPN \
             -i "$DIR/A2D5.idx" \
             -l A \
             -r "$ROOT/${libs[$name]}" \
             -o "$DIR/seed/$name"
