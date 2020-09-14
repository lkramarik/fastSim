#!/bin/bash
mkdir -p workDir
productionId=`date +%F_%H-%M`
echo ${productionId}
mkdir -p workDir/${productionId}
cd workDir/${productionId}
echo "Start of cps."
cp -Lr ../../submit/* ./
cp ../../toyMcEffZeroDecayLength.C ./
cp ../../pp200_spectra.root ./

cp ../../dcaxy_vs_dcaz.root ./
cp ../../tupleManage/divide_ntp.cxx ./

cp ../../totalEff_pi.root ./
cp ../../totalEff_K.root ./

cp ../../kplus_tpc_eff_embedding.root ./
cp ../../kminus_tpc_eff_embedding.root  ./
cp ../../piminus_tpc_eff_embedding.root ./
cp ../../piplus_tpc_eff_embedding.root ./
cp ../../HIJING_D0_pt_y.root ./

cp ../../kaon_momentum_resolution.root ./
cp ../../pion_momentum_resolution.root ./
cp ../../inputs.event.root ./
cp ../../vertexReso.root ./
cp ../../eff_tof.root ./
cp ../../hftratio_vs_pt_dAu_kaon.root ./
cp ../../hftratio_vs_pt_dAu_pion.root ./
cp ../../published_run10_D0_AuAu_data.root ./
cp ../../StAnaCutsData.h ./

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