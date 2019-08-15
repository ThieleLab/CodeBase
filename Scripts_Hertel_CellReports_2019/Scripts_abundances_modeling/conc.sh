########################
#!user/bin/bash
arr=(/media/sf_Y_DRIVE/Microbiome/PD/mapping/PD/ERS*)
# iterate through array using a counter
for ((i=0; i<${#arr[@]}; i++)); do
    cd ${arr[$i]}
    farr1=(./*_1.fastq.gz)
    for ((j=1; j<${#farr1[@]}; j++)) ; do
       cat ${farr1[$j]} >> ${farr1[0]}
    done
    gunzip ${farr1[0]}
    mv ${farr1[0]/".gz"/""} ${arr[$i]##*/}"_1.fastq"

    farr2=(./*_2.fastq.gz)
    for ((j=1; j<${#farr2[@]}; j++)) ; do
       cat ${farr2[$j]} >> ${farr2[0]}
    done
    gunzip ${farr2[0]}
    mv ${farr2[0]/".gz"/""} ${arr[$i]##*/}"_2.fastq"
done
