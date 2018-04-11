#!/bin/bash
mkdir -p workDir
productionId=`date +%F_%H-%M`
echo ${productionId}
mkdir -p workDir/${productionId}
cd workDir/${productionId}
echo "Start of cps."
cp -Lr /gpfs01/star/pwg/lkramarik/sim/fastsimD0_2016/submit/* ./
cp ../../toyMcEffZeroDecayLength.C ./
cp ../../Momentum_resolution_SL16d.root ./
cp ../../pp200_spectra.root ./
cp ../../Run14_D0_MyRaa_pT1.0.root ./
cp ../../HftRatio_AuAu2016_lumiprod.root ./
cp ../../Dca2D_AuAu2016_lumiprod.root ./
cp ../../Vz_Cent.root ./

cp ../../kplus_tpc_eff_embedding.root ./
cp ../../kminus_tpc_eff_embedding.root  ./
cp ../../piminus_tpc_eff_embedding.root ./
cp ../../piplus_tpc_eff_embedding.root ./

cp ../../kaon_momentum_resolution.root ./
cp ../../pion_momentum_resolution.root ./
cp ../../inputs.event.root ./
cp ../../eff_tof.root ./
cp ../../hftratio_vs_pt_dAu_kaon.root ./
cp ../../hftratio_vs_pt_dAu_pion.root ./

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

star-submit-template -template submitToyMcZeroVtx.xml -entities productionId=${productionId},codePath=$codePath



