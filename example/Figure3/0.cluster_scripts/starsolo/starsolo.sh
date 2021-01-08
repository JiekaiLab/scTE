#PBS -N ss.${out}.starsolo
#PBS -l nodes=1:ppn=32
#PBS -l mem=32gb
#PBS -j oe
#PBS -o ss.${out}.out
#PBS -q batch
#PBS -V 
cd $PBS_O_WORKDIR

ulimit -n 2000

whitelist='--soloCBwhitelist /data3/lab-andrew/scTE/scrnaseq_barcodes/version1.txt' # Make sure you get the right bartcode version

# Required arguments;
mods='--soloType Droplet --soloFeatures Gene --soloBarcodeReadLength 1 --soloCBlen 14 --soloUMIstart 15 '
teopts=' --outFilterMultimapNmax 100 --winAnchorMultimapNmax 100 --outSAMmultNmax 1 --outSAMtype BAM SortedByCoordinate --twopassMode Basic' 
opts='--runRNGseed 42 --runThreadN 32 --readFilesCommand zcat '

# required for scTE:
sam_att='--outSAMattributes NH HI AS nM CR CY UR UY'

genome_mm10='--genomeDir /data3/lab-andrew/scTE/custom_indeces/mm10_gencode_vM21_starsolo/SAindex'
genome_hg38='--genomeDir /data3/lab-andrew/scTE/custom_indeces/hg38_gencode_v30_starsolo/SAindex'

# p1 = read
# p2 = barcode and UMI
# Make sure you set the correct genome index;
STAR $opts $teopts $mods $whitelist $sam_att $genome_mm10 --outFileNamePrefix ss.${out} --readFilesIn ${p1} ${p2}

rm -r ss.${out}_STARgenome
rm -r ss.${out}_STARpass1
rm -r ss.${out}_STARtmp
