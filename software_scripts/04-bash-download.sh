
set -x

PATH=$PATH:`pwd`/software/sratoolkit/bin:`pwd`/software/enaBrowserTools/python

cat wheat_98 | while read i
do
enaDataGet -f fastq ${i}
if [ ! -f ${i}.fastq.gz ];then
	enaDataGet -f sra ${i}
if [ ! -f ${i} ];then
	wget ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/${i:0:6}/${i}/${i}.sra
fi
fi
done
