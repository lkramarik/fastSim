#!/bin/bash
mkdir -p workDir
productionId=`date +%F_%H-%M`
echo ${productionId}
mkdir -p workDir/${productionId}
cd workDir/${productionId}
echo "Start of cps."
cp -Lr /gpfs01/star/pwg/lkramarik/sim/fastsimD0_2016/submit/* ./
cp ../../toyMcEffZeroDecayLength.C ./
cp ../../pp200_spectra.root ./
cp ../../Run14_D0_MyRaa_pT1.0.root ./

cp ../../dcaxy_vs_dcaz_hijing.root ./
cp ../../tupleManage/divide_ntp.cxx ./

cp ../../totalEff_pi.root ./
cp ../../totalEff_K.root ./

cp ../../kplus_tpc_eff_embedding.root ./
cp ../../kminus_tpc_eff_embedding.root  ./
cp ../../piminus_tpc_eff_embedding.root ./
cp ../../piplus_tpc_eff_embedding.root ./

cp ../../kaon_momentum_resolution.root ./
cp ../../pion_momentum_resolution.root ./
cp ../../inputs.event.root ./
cp ../../vertexReso.root ./
cp ../../eff_tof.root ./
cp ../../hftratio_vs_pt_dAu_kaon_hijing.root ./
cp ../../hftratio_vs_pt_dAu_pion_hijing.root ./
cp ../../published_run10_D0_AuAu_data.root ./

echo "Copying done."

rm -r LocalLibraries.package
rm LocalLibraries.zip

codePath="/gpfs01/star/pwg/lkramarik/sim/fastsimD0_2016/workDir/"${productionId}
echo codePath

mkdir ./jobs
mkdir ./out/
mkdir ./jobs/log/
mkdir ./jobs/err/
mkdir ./jobs/report/
mkdir ./jobs/list
mkdir ./jobs/csh

echo "Directories created. Let's submit some jobs..."

star-submit-template -template submitToyMcZeroVtxHIJING.xml -entities productionId=${productionId},codePath=$codePath