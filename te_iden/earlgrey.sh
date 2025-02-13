#!/bin/bash
#$-q all.q
#$-V
#$-cwd
#$-pe smp 25
#Change the number of required cores everytime, setaria around 15 and saccharum around 25

module load miniconda3
. "/Storage/progs/miniconda3/etc/profile.d/mamba.sh"
mamba activate /Storage/progs/EarlGrey
echo $PERL5LIB

#Only change the following lines
GENOME_FILE="SofficinarumxspontaneumR570_771_v2.0.fa"
WORK_DIR=/Storage/data2/andreza.cunha/workdir
RESULTS_DIR=EGHybridR570
SPECIES=Shybrid
#Do not change after this line

mkdir -p $WORK_DIR
cp $GENOME_FILE $WORK_DIR
cd $WORK_DIR

#Run earl grey with minimum command options
earlGrey -t 25 -g $GENOME_FILE -s $SPECIES -o $RESULTS_DIR
cd $SGE_O_WORKDIR
mv ${WORK_DIR}/${RESULTS_DIR} .
mv $RESULTS_DIR /Storage/data2/andreza.cunha/results/
rm -rf $WORK_DIR
