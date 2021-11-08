### quality control
time kneaddata --input rawdata/$line'_1'.fq.gz --input rawdata/$line'_2'.fq.gz -t 28 -db Database/hg19/Homo_sapiens 
    --output 01Knead/ --output-prefix $line --trimmomatic-options 'SLIDINGWINDOW:4:20 MINLEN:50 LEADING:3 TRAILING:3'
kneaddata --input 01Knead/$line'_paired_1'.fastq --input 01Knead/$line'_paired_2'.fastq -t 28 -db Database/contaminate_db --output 01Decom/ 
    --output-prefix $line --bypass-trim
### Kraken annotation
Name='Archaea Bacteria fungi viral'
cat SampleID.tsv|while read line
do
    echo $line
    for i in $name
    do
        time kraken2  --db ~/Database/$i --threads 56 --use-names --report krakenS/$line'_'$i'.kreport2' --paired 01Decom/$line'_paired_1'.fastq 01Decom/$line'_paired_2'.fastq --output krakenS/$line'_'$i'.kraken'
        time bracken -d ~/Database/$i -t 56 -l D  -i krakenS/$line'_'$i'.kreport2' -o krakenS/$line'_'$i'.D.braken'
        time bracken -d ~/Database/$i -t 56 -l P  -i krakenS/$line'_'$i'.kreport2' -o krakenS/$line'_'$i'.P.braken'
        time bracken -d ~/Database/$i -t 56 -l C  -i krakenS/$line'_'$i'.kreport2' -o krakenS/$line'_'$i'.C.braken'
        time bracken -d ~/Database/$i -t 56 -l O  -i krakenS/$line'_'$i'.kreport2' -o krakenS/$line'_'$i'.O.braken'
        time bracken -d ~/Database/$i -t 56 -l F  -i krakenS/$line'_'$i'.kreport2' -o krakenS/$line'_'$i'.F.braken'
        time bracken -d ~/Database/$i -t 56 -l G  -i krakenS/$line'_'$i'.kreport2' -o krakenS/$line'_'$i'.G.braken'
        time bracken -d ~/Database/$i -t 56 -l S  -i krakenS/$line'_'$i'.kreport2' -o krakenS/$line'_'$i'.S.braken'
    done
done
#### assembly
megahit -1 $line'_paired_1'.fastq -2 $line'_paired_2'.fastq -m 0.9 -t 28 --presets meta-sensitive --out-prefix $line -o assem/$line'_assem'
### Gene Prediction
prodigal -a Prodigal/$line'_prot.faa' -i assem/$line'_assem'/$line'.contigs.fa' -d Prodigal/$line'_nucl.fa' -o Prodigal/$line'_gene.gff' -f gff -p meta -s Prodigal/$line'_stat'.txt

### Non-redundant gene set
cd-hit-est -i genes_all.fa -o genes_all_unique.fa -aS 0.9 -c 0.95 -G 0 -g 0 -T 28 -M 0

### Gene annotation
emapper.py -m diamond --no_annot --no_file_comments --cpu 28 --data_dir Database/eggnog -i $file -o $file
emapper.py --annotate_hits_table CRC_emapper.emapper.seed_orthologs  --no_file_comments --cpu 28 --data_dir Database/eggnog -o CRC_emapper

### Gene abundance
coverm contig -t 28 -r contig_reference.fa --coupled $fq --bam-file-cache-directory coverm_bamout

