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
        time kraken2 --use-names --use-mpa-style --report-zero-counts --threads 28 --db ~/Database/$i --report Kraken2/$line'_'$i'.kreport2' --paired 01Decom/$line'_paired_1'.fastq 01Decom/$line'_paired_2'.fastq --output Kraken2/$line'_'$i'.kraken'
    done
done
#### assembly
megahit -1 $line'_paired_1'.fastq -2 $line'_paired_2'.fastq -m 0.9 -t 28 --presets meta-sensitive --out-prefix $line -o assem/$line'_assem'
### Gene Prediction
prodigal -a Prodigal/$line'_prot.faa' -i assem/$line'_assem'/$line'.contigs.fa' -d Prodigal/$line'_nucl.fa' -o Prodigal/$line'_gene.gff' -f gff -p meta -s Prodigal/$line'_stat'.txt

### Non-redundant gene set
cd-hit-est -i genes_all.fa -o genes_all_unique.fa -aS 0.9 -c 0.95 -G 0 -g 0 -T 28 -M 0

### Gene annotation
time emapper.py -m diamond --no_annot --no_file_comments --cpu 28 --data_dir Database/eggnog -i $file -o $file
time emapper.py --annotate_hits_table CRC_emapper.emapper.seed_orthologs  --no_file_comments --cpu 28 --data_dir Database/eggnog -o CRC_emapper
