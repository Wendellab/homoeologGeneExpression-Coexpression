#!/bin/bash
#PBS -N salmon_flower -l walltime=48:00:00 -t 0-65
DIR='/work/maa146/EffectiveLength/salmon-edgeR/'
unset libs
declare -A libs
libs['LD7-A2-L1']='raw/flowerTimeDataset/LD7-A2-L1_1.fq.gz
raw/flowerTimeDataset/LD7-A2-L1_2.fq.gz'
libs['LD7-A2-L2']='raw/flowerTimeDataset/LD7-A2-L2_1.fq.gz
raw/flowerTimeDataset/LD7-A2-L2_2.fq.gz'
libs['LD7-A2-L3']='raw/flowerTimeDataset/LD7-A2-L3_1.fq.gz
raw/flowerTimeDataset/LD7-A2-L3_2.fq.gz'
libs['LD7-ADs-L1']='raw/flowerTimeDataset/LD7-ADs-L1_1.fq.gz
raw/flowerTimeDataset/LD7-ADs-L1_2.fq.gz'
libs['LD7-ADs-L2']='raw/flowerTimeDataset/LD7-ADs-L2_1.fq.gz
raw/flowerTimeDataset/LD7-ADs-L2_2.fq.gz'
libs['LD7-ADs-L3']='raw/flowerTimeDataset/LD7-ADs-L3_1.fq.gz
raw/flowerTimeDataset/LD7-ADs-L3_2.fq.gz'
libs['LD7-D5-L1']='raw/flowerTimeDataset/LD7-D5-L1_1.fq.gz
raw/flowerTimeDataset/LD7-D5-L1_2.fq.gz'
libs['LD7-D5-L2']='raw/flowerTimeDataset/LD7-D5-L2_1.fq.gz
raw/flowerTimeDataset/LD7-D5-L2_2.fq.gz'
libs['LD7-D5-L3']='raw/flowerTimeDataset/LD7-D5-L3_1.fq.gz
raw/flowerTimeDataset/LD7-D5-L3_2.fq.gz'
libs['LD9-A2-L1']='raw/flowerTimeDataset/LD9-A2-L1_1.fq.gz
raw/flowerTimeDataset/LD9-A2-L1_2.fq.gz'
libs['LD9-A2-L2']='raw/flowerTimeDataset/LD9-A2-L2_1.fq.gz
raw/flowerTimeDataset/LD9-A2-L2_2.fq.gz'
libs['LD9-A2-L3']='raw/flowerTimeDataset/LD9-A2-L3_1.fq.gz
raw/flowerTimeDataset/LD9-A2-L3_2.fq.gz'
libs['LD9-ADs-L1']='raw/flowerTimeDataset/LD9-ADs-L1_1.fq.gz
raw/flowerTimeDataset/LD9-ADs-L1_2.fq.gz'
libs['LD9-ADs-L2']='raw/flowerTimeDataset/LD9-ADs-L2_1.fq.gz
raw/flowerTimeDataset/LD9-ADs-L2_2.fq.gz'
libs['LD9-ADs-L3']='raw/flowerTimeDataset/LD9-ADs-L3_1.fq.gz
raw/flowerTimeDataset/LD9-ADs-L3_2.fq.gz'
libs['LD9-D5-L1']='raw/flowerTimeDataset/LD9-D5-L1_1.fq.gz
raw/flowerTimeDataset/LD9-D5-L1_2.fq.gz'
libs['LD9-D5-L2']='raw/flowerTimeDataset/LD9-D5-L2_1.fq.gz
raw/flowerTimeDataset/LD9-D5-L2_2.fq.gz'
libs['LD9-D5-L3']='raw/flowerTimeDataset/LD9-D5-L3_1.fq.gz
raw/flowerTimeDataset/LD9-D5-L3_2.fq.gz'
libs['LDM-A2-L1']='raw/flowerTimeDataset/LDM-A2-L1_1.fq.gz
raw/flowerTimeDataset/LDM-A2-L1_2.fq.gz'
libs['LDM-A2-L2']='raw/flowerTimeDataset/LDM-A2-L2_1.fq.gz
raw/flowerTimeDataset/LDM-A2-L2_2.fq.gz'
libs['LDM-A2-L3']='raw/flowerTimeDataset/LDM-A2-L3_1.fq.gz
raw/flowerTimeDataset/LDM-A2-L3_2.fq.gz'
libs['LDM-ADs-L1']='raw/flowerTimeDataset/LDM-ADs-L1_1.fq.gz
raw/flowerTimeDataset/LDM-ADs-L1_2.fq.gz'
libs['LDM-ADs-L2']='raw/flowerTimeDataset/LDM-ADs-L2_1.fq.gz
raw/flowerTimeDataset/LDM-ADs-L2_2.fq.gz'
libs['LDM-ADs-L3']='raw/flowerTimeDataset/LDM-ADs-L3_1.fq.gz
raw/flowerTimeDataset/LDM-ADs-L3_2.fq.gz'
libs['LDM-D5-L1']='raw/flowerTimeDataset/LDM-D5-L1_1.fq.gz
raw/flowerTimeDataset/LDM-D5-L1_2.fq.gz'
libs['LDM-D5-L2']='raw/flowerTimeDataset/LDM-D5-L2_1.fq.gz
raw/flowerTimeDataset/LDM-D5-L2_2.fq.gz'
libs['LDM-D5-L3']='raw/flowerTimeDataset/LDM-D5-L3_1.fq.gz
raw/flowerTimeDataset/LDM-D5-L3_2.fq.gz'
libs['SD5-A2-S1']='raw/flowerTimeDataset/SD5-A2-S1_1.fq.gz
raw/flowerTimeDataset/SD5-A2-S1_2.fq.gz'
libs['SD5-A2-S5']='raw/flowerTimeDataset/SD5-A2-S5_1.fq.gz
raw/flowerTimeDataset/SD5-A2-S5_2.fq.gz'
libs['SD5-ADs-S1']='raw/flowerTimeDataset/SD5-ADs-S1_1.fq.gz
raw/flowerTimeDataset/SD5-ADs-S1_2.fq.gz'
libs['SD5-ADs-S3']='raw/flowerTimeDataset/SD5-ADs-S3_1.fq.gz
raw/flowerTimeDataset/SD5-ADs-S3_2.fq.gz'
libs['SD5-D5-S1']='raw/flowerTimeDataset/SD5-D5-S1_1.fq.gz
raw/flowerTimeDataset/SD5-D5-S1_2.fq.gz'
libs['SD5-D5-S3']='raw/flowerTimeDataset/SD5-D5-S3_1.fq.gz
raw/flowerTimeDataset/SD5-D5-S3_2.fq.gz'
libs['SD7-A2-S1']='raw/flowerTimeDataset/SD7-A2-S1_1.fq.gz
raw/flowerTimeDataset/SD7-A2-S1_2.fq.gz'
libs['SD7-A2-S4']='raw/flowerTimeDataset/SD7-A2-S4_1.fq.gz
raw/flowerTimeDataset/SD7-A2-S4_2.fq.gz'
libs['SD7-A2-S5']='raw/flowerTimeDataset/SD7-A2-S5_1.fq.gz
raw/flowerTimeDataset/SD7-A2-S5_2.fq.gz'
libs['SD7-ADs-S1']='raw/flowerTimeDataset/SD7-ADs-S1_1.fq.gz
raw/flowerTimeDataset/SD7-ADs-S1_2.fq.gz'
libs['SD7-ADs-S2']='raw/flowerTimeDataset/SD7-ADs-S2_1.fq.gz
raw/flowerTimeDataset/SD7-ADs-S2_2.fq.gz'
libs['SD7-ADs-S3']='raw/flowerTimeDataset/SD7-ADs-S3_1.fq.gz
raw/flowerTimeDataset/SD7-ADs-S3_2.fq.gz'
libs['SD7-D5-S1']='raw/flowerTimeDataset/SD7-D5-S1_1.fq.gz
raw/flowerTimeDataset/SD7-D5-S1_2.fq.gz'
libs['SD7-D5-S2']='raw/flowerTimeDataset/SD7-D5-S2_1.fq.gz
raw/flowerTimeDataset/SD7-D5-S2_2.fq.gz'
libs['SD7-D5-S3']='raw/flowerTimeDataset/SD7-D5-S3_1.fq.gz
raw/flowerTimeDataset/SD7-D5-S3_2.fq.gz'
libs['SDM-A2-S5']='raw/flowerTimeDataset/SDM-A2-S5_1.fq.gz
raw/flowerTimeDataset/SDM-A2-S5_2.fq.gz'
libs['SDM-A2-S6']='raw/flowerTimeDataset/SDM-A2-S6_1.fq.gz
raw/flowerTimeDataset/SDM-A2-S6_2.fq.gz'
libs['SDM-ADs-S1']='raw/flowerTimeDataset/SDM-ADs-S1_1.fq.gz
raw/flowerTimeDataset/SDM-ADs-S1_2.fq.gz'
libs['SDM-ADs-S2']='raw/flowerTimeDataset/SDM-ADs-S2_1.fq.gz
raw/flowerTimeDataset/SDM-ADs-S2_2.fq.gz'
libs['SDM-D5-S1']='raw/flowerTimeDataset/SDM-D5-S1_1.fq.gz
raw/flowerTimeDataset/SDM-D5-S1_2.fq.gz'
libs['SDM-D5-S2']='raw/flowerTimeDataset/SDM-D5-S2_1.fq.gz
raw/flowerTimeDataset/SDM-D5-S2_2.fq.gz'
libs['SDP9-A2-S1']='raw/flowerTimeDataset/SDP9-A2-S1_1.fq.gz
raw/flowerTimeDataset/SDP9-A2-S1_2.fq.gz'
libs['SDP9-A2-S3']='raw/flowerTimeDataset/SDP9-A2-S3_1.fq.gz
raw/flowerTimeDataset/SDP9-A2-S3_2.fq.gz'
libs['SDP9-A2-S4']='raw/flowerTimeDataset/SDP9-A2-S4_1.fq.gz
raw/flowerTimeDataset/SDP9-A2-S4_2.fq.gz'
libs['SDP9-ADs-S1']='raw/flowerTimeDataset/SDP9-ADs-S1_1.fq.gz
raw/flowerTimeDataset/SDP9-ADs-S1_2.fq.gz'
libs['SDP9-ADs-S2']='raw/flowerTimeDataset/SDP9-ADs-S2_1.fq.gz
raw/flowerTimeDataset/SDP9-ADs-S2_2.fq.gz'
libs['SDP9-ADs-S3']='raw/flowerTimeDataset/SDP9-ADs-S3_1.fq.gz
raw/flowerTimeDataset/SDP9-ADs-S3_2.fq.gz'
libs['SDP9-D5-S1']='raw/flowerTimeDataset/SDP9-D5-S1_1.fq.gz
raw/flowerTimeDataset/SDP9-D5-S1_2.fq.gz'
libs['SDP9-D5-S2']='raw/flowerTimeDataset/SDP9-D5-S2_1.fq.gz
raw/flowerTimeDataset/SDP9-D5-S2_2.fq.gz'
libs['SDP9-D5-S3']='raw/flowerTimeDataset/SDP9-D5-S3_1.fq.gz
raw/flowerTimeDataset/SDP9-D5-S3_2.fq.gz'
libs['SDPF-A2-S1']='raw/flowerTimeDataset/SDPF-A2-S1_1.fq.gz
raw/flowerTimeDataset/SDPF-A2-S1_2.fq.gz'
libs['SDPF-A2-S3']='raw/flowerTimeDataset/SDPF-A2-S3_1.fq.gz
raw/flowerTimeDataset/SDPF-A2-S3_2.fq.gz'
libs['SDPF-A2-S6']='raw/flowerTimeDataset/SDPF-A2-S6_1.fq.gz
raw/flowerTimeDataset/SDPF-A2-S6_2.fq.gz'
libs['SDPF-ADs-S1']='raw/flowerTimeDataset/SDPF-ADs-S1_1.fq.gz
raw/flowerTimeDataset/SDPF-ADs-S1_2.fq.gz'
libs['SDPF-ADs-S2']='raw/flowerTimeDataset/SDPF-ADs-S2_1.fq.gz
raw/flowerTimeDataset/SDPF-ADs-S2_2.fq.gz'
libs['SDPF-ADs-S3']='raw/flowerTimeDataset/SDPF-ADs-S3_1.fq.gz
raw/flowerTimeDataset/SDPF-ADs-S3_2.fq.gz'
libs['SDPF-D5-S1']='raw/flowerTimeDataset/SDPF-D5-S1_1.fq.gz
raw/flowerTimeDataset/SDPF-D5-S1_2.fq.gz'
libs['SDPF-D5-S2']='raw/flowerTimeDataset/SDPF-D5-S2_1.fq.gz
raw/flowerTimeDataset/SDPF-D5-S2_2.fq.gz'
libs['SDPF-D5-S3']='raw/flowerTimeDataset/SDPF-D5-S3_1.fq.gz
raw/flowerTimeDataset/SDPF-D5-S3_2.fq.gz'
cd $DIR/flower
ROOT=$(git rev-parse --show-toplevel)

PATH=$PATH:$DIR/Salmon-latest_linux_x86_64/bin

KEYS=("${!libs[@]}")
name="${KEYS[$PBS_ARRAYID]}"

readarray -t files <<< "${libs[$name]}"

salmon quant --no-version-check \
             -p $PBS_NUM_PPN \
             -i "$DIR/A2D5.idx" \
             -l A \
             -1 "$ROOT/${files[0]}" \
             -2 "$ROOT/${files[1]}" \
             -o "$DIR/flower/$name"
