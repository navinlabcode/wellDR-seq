WPATH=`pwd`
RNAPATH=$WPATH/rawdata
REF_Genome=/volumes/USR2/wangkl/db/star_2.7.5_hg19_3.0.0
ANNO=/volumes/seq/code/3rd_party/10X/refdata-cellranger-hg19-3.0.0/genes/genes.gtf
DATE=$(date +"%Y-%m-%d")
OUT_PATH=FILE_exact_hg19_star_smart_2
NUM_thread=40
WHITELIST=/volumes/USR2/wangkl/wscDR/barcode/com_index_wafer.txt
RAWPATH=$(echo $RNAPATH | sed 's/USR2\/wangkl/USR\/INACTIVE_users\/kaile_Data/')
TRIM_DIR=trimmed_FILE_R2_fq_list_rna_2
source ~/.bashrc
source ~/.bash_profile

cd $RAWPATH
#:<<B1
start_time=`date +%s`
echo "FILE scDR-fq-convert start at: `date`"
perl ~/bin/scDR-fq-convert-v4.pl FILE_R1.fastq.gz FILE_R2.fastq.gz $WHITELIST $NUM_thread
echo "FILE scDR-fq-convert finish at: `date`"

mkdir $TRIM_DIR
echo "FILE Trim stat at: `date`"
#---Multi-thread Trim
split -l $((5184/$NUM_thread)) $WHITELIST FILE_whlst.
function trim_bg () {
for i in $(cat $j)
do
java -jar ~/bin/Trimmomatic-0.39/trimmomatic-0.39.jar SE -threads $NUM_thread ./FILE_R2_fq_list_rna_m/$i\.fq.gz ./$TRIM_DIR/trimmed_$i\.fq.gz ILLUMINACLIP:/volumes/USR2/wangkl/wscDR/trim_adp.fa:2:30:10 LEADING:0 TRAILING:0 SLIDINGWINDOW:4:0 MINLEN:20 >> FILE_smart_trim.log 2>&1
done
}
for j in $(ls FILE_whlst.*)
do
trim_bg $j &
done
wait
rm FILE_whlst.*
echo "FILE Trim finish at: `date`"
#B1
#:<<B1
awk -v rp="$RAWPATH" -v trm="$TRIM_DIR" '{print rp"/"trm"/trimmed_"$1".fq.gz\t-\t"$1}' $WHITELIST > star_manifest_FILE.tsv

cd $WPATH/data/
mkdir -p $OUT_PATH/star_2.7.5
#B1
cd $WPATH/data/$OUT_PATH/star_2.7.5
#:<<B1
/volumes/seq/code/3rd_party/STAR/STAR-2.7.5a/bin/Linux_x86_64/STAR --runThreadN $NUM_thread --genomeDir $REF_Genome --readFilesCommand zcat --outFilterMultimapNmax 1 --outSAMtype BAM SortedByCoordinate --outFileNamePrefix $OUT_PATH\_ --soloType SmartSeq --readFilesManifest $RAWPATH/star_manifest_FILE.tsv --soloUMIdedup Exact NoDedup --soloStrand Unstranded --soloFeatures Gene GeneFull

mkdir ./$OUT_PATH\_Solo.out/Gene/filtered_exact/
mkdir ./$OUT_PATH\_Solo.out/Gene/filtered_nodup/
mkdir ./$OUT_PATH\_Solo.out/GeneFull/filtered_exact/
mkdir ./$OUT_PATH\_Solo.out/GeneFull/filtered_nodup/
mv ./$OUT_PATH\_Solo.out/Gene/filtered/matrix.mtx ./$OUT_PATH\_Solo.out/Gene/filtered/matrix.mtx_old
mv ./$OUT_PATH\_Solo.out/GeneFull/filtered/matrix.mtx ./$OUT_PATH\_Solo.out/GeneFull/filtered/matrix.mtx_old
#B1
cp ./$OUT_PATH\_Solo.out/Gene/filtered/* ./$OUT_PATH\_Solo.out/Gene/filtered_exact/
cp ./$OUT_PATH\_Solo.out/Gene/filtered/* ./$OUT_PATH\_Solo.out/Gene/filtered_nodup/
cp ./$OUT_PATH\_Solo.out/GeneFull/filtered/* ./$OUT_PATH\_Solo.out/GeneFull/filtered_exact/
cp ./$OUT_PATH\_Solo.out/GeneFull/filtered/* ./$OUT_PATH\_Solo.out/GeneFull/filtered_nodup/
#B1
head -n 3 ./$OUT_PATH\_Solo.out/Gene/filtered/matrix.mtx_old > ./$OUT_PATH\_Solo.out/Gene/filtered/matrix.head
awk 'NR>3' ./$OUT_PATH\_Solo.out/Gene/filtered_exact/matrix.mtx_old | awk '{print $1" "$2" "$4}' | cat ./$OUT_PATH\_Solo.out/Gene/filtered/matrix.head - > ./$OUT_PATH\_Solo.out/Gene/filtered_exact/matrix.mtx
awk 'NR>3' ./$OUT_PATH\_Solo.out/Gene/filtered_nodup/matrix.mtx_old | awk '{print $1" "$2" "$3}' | cat ./$OUT_PATH\_Solo.out/Gene/filtered/matrix.head - > ./$OUT_PATH\_Solo.out/Gene/filtered_nodup/matrix.mtx

head -n 3 ./$OUT_PATH\_Solo.out/GeneFull/filtered/matrix.mtx_old > ./$OUT_PATH\_Solo.out/GeneFull/filtered/matrix.head
awk 'NR>3' ./$OUT_PATH\_Solo.out/GeneFull/filtered_exact/matrix.mtx_old | awk '{print $1" "$2" "$4}' | cat ./$OUT_PATH\_Solo.out/GeneFull/filtered/matrix.head - > ./$OUT_PATH\_Solo.out/GeneFull/filtered_exact/matrix.mtx
awk 'NR>3' ./$OUT_PATH\_Solo.out/GeneFull/filtered_nodup/matrix.mtx_old | awk '{print $1" "$2" "$3}' | cat ./$OUT_PATH\_Solo.out/GeneFull/filtered/matrix.head - > ./$OUT_PATH\_Solo.out/GeneFull/filtered_nodup/matrix.mtx

gzip ./$OUT_PATH\_Solo.out/Gene/filtered_exact/*
gzip ./$OUT_PATH\_Solo.out/Gene/filtered_nodup/*
gzip ./$OUT_PATH\_Solo.out/GeneFull/filtered_exact/*
gzip ./$OUT_PATH\_Solo.out/GeneFull/filtered_nodup/*

mkdir ./$OUT_PATH\_Solo.out/Gene/raw_exact/                                                                                                                                                                                            
mkdir ./$OUT_PATH\_Solo.out/Gene/raw_nodup/                                                                                                                                                                                            
mkdir ./$OUT_PATH\_Solo.out/GeneFull/raw_exact/                                                                                                                                                                                        
mkdir ./$OUT_PATH\_Solo.out/GeneFull/raw_nodup/                                                                                                                                                                                        
mv ./$OUT_PATH\_Solo.out/Gene/raw/matrix.mtx ./$OUT_PATH\_Solo.out/Gene/raw/matrix.mtx_old                                                                                                                                        
mv ./$OUT_PATH\_Solo.out/GeneFull/raw/matrix.mtx ./$OUT_PATH\_Solo.out/GeneFull/raw/matrix.mtx_old                                                                                                                                
#B1                                                                                                                                                                                                                                         
cp ./$OUT_PATH\_Solo.out/Gene/raw/* ./$OUT_PATH\_Solo.out/Gene/raw_exact/                                                                                                                                                         
cp ./$OUT_PATH\_Solo.out/Gene/raw/* ./$OUT_PATH\_Solo.out/Gene/raw_nodup/                                                                                                                                                         
cp ./$OUT_PATH\_Solo.out/GeneFull/raw/* ./$OUT_PATH\_Solo.out/GeneFull/raw_exact/                                                                                                                                                 
cp ./$OUT_PATH\_Solo.out/GeneFull/raw/* ./$OUT_PATH\_Solo.out/GeneFull/raw_nodup/                                                                                                                                                 
#B1
head -n 3 ./$OUT_PATH\_Solo.out/Gene/raw/matrix.mtx_old > ./$OUT_PATH\_Solo.out/Gene/raw/matrix.head
awk 'NR>3' ./$OUT_PATH\_Solo.out/Gene/raw_exact/matrix.mtx_old | awk '{print $1" "$2" "$4}' | cat ./$OUT_PATH\_Solo.out/Gene/raw/matrix.head - > ./$OUT_PATH\_Solo.out/Gene/raw_exact/matrix.mtx
awk 'NR>3' ./$OUT_PATH\_Solo.out/Gene/raw_nodup/matrix.mtx_old | awk '{print $1" "$2" "$3}' | cat ./$OUT_PATH\_Solo.out/Gene/raw/matrix.head - > ./$OUT_PATH\_Solo.out/Gene/raw_nodup/matrix.mtx

head -n 3 ./$OUT_PATH\_Solo.out/GeneFull/raw/matrix.mtx_old > ./$OUT_PATH\_Solo.out/GeneFull/raw/matrix.head
awk 'NR>3' ./$OUT_PATH\_Solo.out/GeneFull/raw_exact/matrix.mtx_old | awk '{print $1" "$2" "$4}' | cat ./$OUT_PATH\_Solo.out/GeneFull/raw/matrix.head - > ./$OUT_PATH\_Solo.out/GeneFull/raw_exact/matrix.mtx
awk 'NR>3' ./$OUT_PATH\_Solo.out/GeneFull/raw_nodup/matrix.mtx_old | awk '{print $1" "$2" "$3}' | cat ./$OUT_PATH\_Solo.out/GeneFull/raw/matrix.head - > ./$OUT_PATH\_Solo.out/GeneFull/raw_nodup/matrix.mtx

gzip ./$OUT_PATH\_Solo.out/Gene/raw_exact/*
gzip ./$OUT_PATH\_Solo.out/Gene/raw_nodup/*
gzip ./$OUT_PATH\_Solo.out/GeneFull/raw_exact/*
gzip ./$OUT_PATH\_Solo.out/GeneFull/raw_nodup/*

chmod -R 755 $WPATH/data
