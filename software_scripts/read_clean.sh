set -x
# phredval=`perl software_scripts/fastq_phred.pl $1 | grep -o "Phred.." | sed "s/Phred//" `

# yn_barcode=`head -n 1000000 $1 | \
# awk '{if((NR+2)%4==0){a[substr($1,1,2)]++;b[substr($1,1,3)]++;c[substr($1,1,4)]++}} \
# END{if(length(c)==1){print 5}else if(length(b)==1){print 4} \
# else if(length(a)==1){print 3}else{print "No"}}' `
# yn_adapter=`head -n 1000000 $1 | grep "length=21"  | wc -l`
# if [ $yn_barcode != "No" ];then
#     software_scripts/fastx/fastx_trimmer -f ${yn_barcode} -i $1 -o $2
#     mv $2 $1
# fi
# if [ $yn_adapter -lt 1 ];then
#     adapter=`python software_scripts/0-dnapi.py $1`
#     software_scripts/fastx/fastx_clipper -n -v -Q${phredval} -l 10 -a $adapter -i $1 -o $2
#     mv $2 $1
# fi
# software_scripts/fastx/fastq_quality_filter -q 20 -p 80 -Q${phredval} -i ${1}.tmp -o $2
# echo "$3\t${yn_barcode}\t${adapter}\t${phredval}" >>$4

adapter=`grep $3 Adapter.txt | cut -f 2`
software_scripts/fastx/fastx_clipper -n -v -Q33 -l 10 -a $adapter -i $1 -o ${1}.tmp1
software_scripts/fastx/fastq_quality_filter -q 20 -p 80 -Q33 -i ${1}.tmp1 -o ${1}.tmp2
software_scripts/fastx/fastq_to_fasta -r -n -v -Q33 -i ${1}.tmp2 -o $2
rm ${1}.tmp1 ${1}.tmp2
