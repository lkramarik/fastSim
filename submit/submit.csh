#!/bin/csh


rm -r LocalLibraries.package
rm LocalLibraries.zip

set prodId=`date +%F_%H-%M`
codePath="/gpfs01/star/pwg/lkramarik/sim/fastsimD0_2016/workDir/"

mkdir ./myOutput/${prodId}
mkdir ./jobs/log/${prodId}
mkdir ./jobs/err/${prodId}


star-submit-template -template submitToyMcZeroVtx.xml -entities productionId=${prodId},codePath=$codePath