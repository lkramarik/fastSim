#!/bin/sh
mxEve=50000
codePath=/global/project/projectdirs/starprod/rnc/rsooraj/fastsimuD0_2016/runfastsimD0_2016
outPath=/global/project/projectdirs/starprod/rnc/rsooraj/fastsimuD0_2016/runfastsimD0_2016
inpPath=/star/data36/reco/AuAu_200_production_2016/ReversedFullField/P16ij/2016/
maxFiles=30 #defines number of Wtree files per job (specify indep. by sample)

sample=run16_qa
echo outPath=$outPath

star-submit-template -p pdsf_slurm_chos -template ./multiTreeReaderJobTempl.xml -entities  n1=$mxEve,outPath=$outPath,inpPath=$inpPath,codePath=$codePath,sample=$sample,maxFiles=$maxFiles

 
