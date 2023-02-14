#!/bin/bash/

if [ ! -d input/genotype/plink ]; then mkdir input/genotype/plink; fi
if [ ! -d pipeline/1.4.kinship ]; then mkdir pipeline/1.4.kinship; fi

rsync -xqr /projects/PPC/analysis/ppc_eqtls/pipeline/1.3.genotype/pca/plink/merged.* input/genotype/plink/.

scripts/plink --bfile input/genotype/plink/merged --make-rel square --out pipeline/1.4.kinship/kinship
