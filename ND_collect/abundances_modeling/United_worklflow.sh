

#################----------------------connectiong with MeFit
cd /
cd home/federico
cd lsvirtualenv/local/bin
source ./activate
cd casper_v0.8.2
#./mefit -s sample -r1 sample_R1.fastq -r2 sample_R2.fastq -avgq 


##Now defining input variables
seqPath="/media/sf_Y_DRIVE/Federico/NDcollect/Compressed"
outDir="/media/sf_Y_DRIVE/Federico/NDcollect/Merged" 
MeFitPath="/home/federico/lsvirtualenv/bin/casper_v0.8.2"
run=$MeFitPath"/mefit"  
nwok=3

#To subselect elements of list
cd $seqPath
dirlist=('*R1_001.fastq')
set -- $dirlist

#Analysis and output in right folder
for i in $dirlist ; do
name=$(echo $i | cut -f1 -d_)_$(echo $i | cut -f2 -d_)
#echo $name
outnam=$name".fastq"
output=$outDir"/"$outnam
inputO=$seqPath"/"$name"_L001_R1_001.fastq"
inputT=$seqPath"/"$name"_L001_R2_001.fastq"
cd $outDir
#./mefit -s $output -r1 $inputO -r2 $inputT -avgq 20
$run -s $name -r1 $inputO -r2 $inputT -avgq 20
cd $seqPath
done

#################----------------------connectiong with SPINGO
##converting and Classifiying (SPINGO)
outDir="/media/sf_Y_DRIVE/Federico/NDcollect/Classified"
spingPath="/home/federico/SPINGO/spingo"
mergedDir="/media/sf_Y_DRIVE/Federico/NDcollect/Merged"
nwok=3

#To subselect elements of list
cd $seqPath
dirlist=('*R1_001.fastq')
set -- $dirlist

#classifier and output in right folder
for i in $dirlist ; do
name=$(echo $i | cut -f1 -d_)_$(echo $i | cut -f2 -d_)
echo $name
outnam=$name".fsa"
#input=$mergedDir"/"$name".fastq"
input=$mergedDir"/"$name".ovlp.hq.fastq"
output=$outDir"/"$outnam
sed '/^@/!d;s//>/;N' $input> $output #converting in fasta
input=$output
outnam=$name".out"
output=$outDir"/"$outnam
~/SPINGO/spingo -p $nwok -d ~/SPINGO/database/RDP_11.2.species.fa -i $input > $output
done
