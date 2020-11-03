#!/bin/bash
# Usage:  makedagfile.sh <infile> > mydagfile.dag
# Create processing dagfile for each filebase line in file <infile>
# 
if [ $# -ne 1 ]; then
    echo "Incorrect number of arguments."
    echo "Usage:  makedagfile.sh <infile> > mydagfile.dag"    
    exit 1
fi


while IFS= read -r line ; do
   filebase="$line"
echo '#'  $filebase
echo job sam-$filebase submit-1sam.jdl 
echo job bed-$filebase submit-2bed.jdl
echo job star-$filebase submit-3star.jdl
echo job sort-$filebase submit-4sort.jdl
echo job gatk-$filebase  submit-5gatk.jdl
echo job clean-$filebase  submit-6clean.jdl
echo vars sam-$filebase filebase=\"$filebase\"
echo vars bed-$filebase filebase=\"$filebase\"
echo vars star-$filebase filebase=\"$filebase\"
echo vars sort-$filebase filebase=\"$filebase\"
echo vars gatk-$filebase filebase=\"$filebase\"
echo vars clean-$filebase filebase=\"$filebase\"
echo parent sam-$filebase child bed-$filebase
echo parent bed-$filebase child star-$filebase
echo parent star-$filebase child sort-$filebase
echo parent sort-$filebase child gatk-$filebase
echo parent gatk-$filebase child clean-$filebase
echo ""
done < $1

echo vars all_nodes indir=\"/data/GTEx/v8_jon/bams\"
echo vars all_nodes genomedir=\"/data/genomes/GRCh38_Gencode25\"
echo vars all_nodes setup=\"/data/jwerner/tools/setup.sh\"
echo vars all_nodes outdir=\"/data/hover/work/werner1\"


