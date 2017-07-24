# Thesis_Workflow

Begin by downloading and unpacking the genomic data from www.Xenbase.org

The Xenopus laevis JGI 9.1 genome
```
wget ftp://ftp.xenbase.org/pub/Genomics/JGI/Xenla9.1/Xla.v91.repeatMasked.fa.gz
  
gunzip Xla.v91.repeatMasked.fa.gz
```
The X.laevis 9.1 primary transcripts for use in the Salmon pipeline
```
wget ftp://ftp.xenbase.org/pub/Genomics/JGI/Xenla9.1/1.8.3.2/XL_9.1_v1.8.3.2.primaryTranscripts.fa.gz 

gunzip XL_9.1_v1.8.3.2.primaryTranscripts.fa.gz 
```
The X.laevis gff3 file. Create a GTF file using this.
```
wget ftp://ftp.xenbase.org/pub/Genomics/JGI/Xenla9.1/1.8.3.2/XL_9.1_v1.8.3.2.primaryTranscripts.gff3.gz

gunzip XL_9.1_v1.8.3.2.primaryTranscripts.gff3.gz

gffread XL_9.1_v1.8.3.2.primaryTranscripts.gff3 -T -o my.gtf
```
Next step is to download and unpack the RNA-seq data from the SRA
```
for i in {51..66};
do
wget ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR3370/SRR33740"$i"/SRR33740"$i".sra -P /output/directory
done

for i in {51..66};
do
fastq-dump --split-files /path/to/.sra/files/SRR33740"$i" -o /path/to/output/directory/ 
done
```
Quality inspection of the raw data using fastqc and MultiQC.
```
files=/path/to/raw/fastq/files/*.fastq

for i in $files;
do 
fastqc -f fastq -o /path/to/output/directory $i
done

conda install -c bioconda multiqc
cd /path/to/fastqc/.zip/files
multiqc ./
firefox multiqc_report.html
```
TrimGalore quality control and inspection with MultiQC
```
dir='/path/to/unprocessed/fastq/files'
for i in {51..66};
do
/path/to/trim_galore -q 30 -length 36 --paired --fastqc \
-o /path/to/output/directory \
$dir/SRR33740"$i"_1.fastq \
$dir/SRR33740"$i"_2.fastq
done

cd /path/to/fastqc/.zip/files
multiqc ./
firefox multiqc_report.html
```
Salmon Indexing and quanitification.
```
salmon index -t /path/to/treanscript/file/transcripts.fa \
-i /path/to/output/directory/xl_index_salmon \
-k 31

for i in {51..66}
do
salmon quant -i /path/to/salmon/index/xl_index_salmon -l A \
-1 /path/to/first/processed/fastq/read/SRR33740"$i"_1_val_1.fq \
-2 /path/to/second/processed/fastq/read/SRR33740"$i"_2_val_2.fq \
-o /path/to/output/directory/SRR33740"$i" \
--seqBias -p 6
done
```
Once quantification has completed the output is modified so that it conforms with the transcript to gene file.
This modification is done using the quants.py script
