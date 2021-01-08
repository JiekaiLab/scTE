#PBS -l nodes=1:ppn=2,mem=64gb
#PBS -j oe
#PBS -o ${out}.out
#PBS -q batch
#PBS -V 
cd $PBS_O_WORKDIR

genome_mm10='/data3/lab-andrew/scTE/scte_indeces/mm10.exclusive.idx'
genome_hg38='/data3/lab-andrew/scTE/scte_indeces/hg38.exclusive.idx'

python3 /share/apps/genomics/unstable/scTE/bin/scTE -i ${in} -x $genome_mm10 -g mm10  -p 1 -o ${out}

gzip ${out}.csv
