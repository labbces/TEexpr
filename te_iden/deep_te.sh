#!/bin/bash

#First script DeepTE_domain helps to find LTR domains and not missclassify
for file in TEs_grasses/*.flTE.fa; do
	outdir_domain=${file/.flTE.fa/_deepTEdomain}
	outdir_results=${file/.flTE.fa/_deepTEresults}
	mkdir -p $outdir_domain
	mkdir -p $outdir_results
	python3 DeepTE_domain.py -d $outdir_domain -o $outdir_domain -i $file -s supfile_dir --hmmscan /home/user/miniconda3/envs/py36/bin/hmmscan
#Second script Deep_TE.py runs final results considering DeepTE_comain.py output
	python3 DeepTE.py -d $outdir_results -o $outdir_results -i $file -sp P -m_dir Plants_model/ -modify $outdir_domain/opt_te_domain_pattern.txt
done
