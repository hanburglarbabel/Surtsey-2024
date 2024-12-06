# 16S rRNA sequencing data analysis

## For non-original data:
(skip to pipeline for local sequencing data)

### download data from NCBI/ENA archives
```
fasterq-dump SRRXXXXXXXXXX SRRXXXXXXXX ....
```
##### unzip fastq.gz files
```
gunzip *.gz	
```
---
## FastQC
- open fastq files in FastQC to check quality and at what length quality deteriorates to determine the length to trim. (if using downloaded data check for adapter content and add to cutadapt if additional adapters need removed)

## Setup 
- create directory to work in, name this folder after your data
```bash
mkdir XX
cd XX
```
- move fastq files, **silvamod138.fasta**, **filter_seq_by_OTUtable.py**, and **LULU.R** to the working directory

From raw sequences pulled from NELS CGB5 and blank from CGB6(tag 70)
converted bam to fastq
```bash
ls *bam | parallel -j 2 "samtools fastq {} > {.}.fastq"
```
# Pipeline
VSEARCH pipeline reference
https://github.com/torognes/vsearch/wiki/Alternative-VSEARCH-pipeline
https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5075697/
_______________________________________________________________________
Open python terminal to run script
________________________________________________________________________
# START OF PIPELINE
``` bash
conda activate pipeline16s

! /bin/bash

VSEARCH=$(which vsearch)
THREADS=6

---
# Amount reads in fastq file.  
for f in *.fastq; do  
 awk '{s++}END{print s/4}' $f  
done
```
### Remove primers 
Change forward primer to the primer sequence used. 

### Trim forward primers
Jørgensen and Zhao: V4
F: 519f (5′-CAGCMGCCGCGGTAA-3')
R: 805r (5′-GACTACHVGGGTATCTAATCC-3')
Home: (Babel, Turke, ext) V4
F: 519f (5'-CAGCMGCCGCGGTAA-3')
R: 805r (5’- GACTACHVGGGTATCTAATCC-3’)

Bergsten: V4
F: 515f (5'-GTGCCAGCMGCCGCGGTAA-3')
R: 806r (5'-GGACTACHVGGGTWTCTAAT-3')

Lee: V4
F: 515f (5′-GTGCCAGCMGCCGCGGTAA-3′)
R: 806r (5′-GGACTACHVGGGTWTCTAAT-3′)

Motamedi: V4 
F: 515f (5'-GTGCCAGCMGCCGCGGTAA-3')
R: 806r (5'-TAATCTWTGGGVHCATCAGG-3')

Jones/Orcutt: V4-V5
F: 515f (5'-GTGYCAGCMGCCGCGGTAA-3')
R: 926r (5'-CCGYCAATTYMTTTRAGTTT-3')

Wee: V4-V5
F: 515f (5'-GTGYCAGCMGCCGCGGTAA-3')
R: 926R (5′‐CCGYCAATTYMTTTRAGTTT-3')

Ramirez: V1-V3
F: 28F (5′- GAGTTTGATCNTGGCTCA G -3′)
R: 388R (5’- TGCTGCCTCCCGTAGGAGT -3’)

Change -g to appropriate forward primer sequence
``` bash
for f in *.fastq; do
    s=$(awk 'BEGIN{FS=OFS="."}{NF--; print}' <<< $f)
     echo
     echo ====================================
     echo Processing sample $s
     echo ====================================    
     s=$(awk -F'.' '{print "2-" $1}' <<< $s)
     echo $s

     # Primers used.
     # Forward primer: 519f (5-́###-3)
     # Referse primer: 805r (5-́###-3)


 cutadapt -j 3 -g CAGCMGCCGCGGTAA --no-indels --discard-untrimmed -o $s.noprimers.fastq $f
 
 done
```

Motamedi had something called dual-indexed primers: 
#```bash
#cutadapt --discard-trimmed --error-rate 0.10 -g ^TTAGAWACCCVHGTAGTCCGGCTGACTGACT -o $s.noprimers.fastq $f done
#```

Print number of reads:
```bash
# Amount reads in fastq file.  
for f in 2-*.fastq; do  
 awk '{s++}END{print s/4}' $f  
done
```
---

### Trim the sequences at 220bp.  
 This can be shortened or lengthened based on what you saw in FastQC. 220 is a good length usually. 
``` bash
 
 for f in 2-*.fastq; do
     s=$(awk 'BEGIN{FS=OFS="."}{NF--; print}' <<< $f)
     echo
     echo ====================================
     echo Processing sample $s
     echo ====================================    
     s=$(awk -F'-' '{print "3-" $2}' <<< $s)

     $VSEARCH --threads $THREADS \
         --fastq_filter $f \
         --fastq_maxns 0 \
         --fastq_trunclen 220 \
         --fastqout $s.trimmed.fastq
 done

# Amount reads in fastq file.  
for f in 3-*.fastq; do  
 awk '{s++}END{print s/4}' $f  
done

```
---
 
### Quality filtering at maxee = 2.
 This sets the maximum errors allowed. It correlates the the Qulaity/Phred score given to each base - this score is the probability the base is incorrect. 
 https://www.drive5.com/usearch/manual/exp_errs.html 
default is 2, but this number can be lowered for more aggresive filtering 

``` bash
  for f in 3-*.fastq; do
     s=$(awk 'BEGIN{FS=OFS="."}{NF--; print}' <<< $f)
     echo
     echo ====================================
     echo Processing sample $s
     echo ====================================    
     s=$(awk -F'-' '{print "4-" $2}' <<< $s)
     l=$(awk -F'-' '{print $2}' <<< $s)
     l=$(awk -F'.' '{print $1}' <<< $l)
     $VSEARCH --threads $THREADS \
         --fastq_filter $f \
         --relabel $l. \
         --fastq_maxee 2 \
         --fastaout $s.fasta \
         --fasta_width 0
 done

for f in 4-*.fasta; do  
 awk '{s++}END{print s/2}' $f  
done
```
---

### Dereplicate at sample level and relabel with sample_n
This removes identical sequences within a sample
``` bash
 for f in 4-*.fasta; do
     s=$(awk 'BEGIN{FS=OFS="."}{NF--; print}' <<< $f)
     echo
     echo ====================================
     echo Processing sample $s
     echo ====================================    
     s=$(awk -F'-' '{print "5-" $2}' <<< $s)
     l=$(awk -F'-' '{print $2}' <<< $s)
     l=$(awk -F'.' '{print $1}' <<< $l)
     $VSEARCH --threads $THREADS \
         --derep_fulllength $f \
         --strand plus \
         --sizeout \
         --relabel $l. \
         --fasta_width 0 \
         --output $s.derep.dfa
 done 

for f in 5-*.dfa; do  
 awk '{s++}END{print s/2}' $f  
done
```
---
### Cat all fasta files
combine all files into one fasta file
``` bash
cat 5-*.dfa > all.fasta
```
---
### Remove unneeded files.
``` bash
rm 2-*.fastq 3-*.fastq 4-*.fasta 5-*.dfa
```
---
### Dereplication of the fasta file.
 removes identical sequences 

**derep.uc** is a way of mapping which sequences goes to which sample and how many times
 ``` bash 
 $VSEARCH --derep_fulllength all.fasta \
     --threads $THREADS \
     --minuniquesize 2 \
     --sizein \
     --sizeout \
     --fasta_width 0 \
    --uc derep.uc \
    --output derep.fasta
```
---
### Clustering at 97% similarity.
De Novo clustering of  sequences with 97% similarity as default (similarity level can be adjusted).  See Rognes et. al (2016) for more in depth explanation
 ```bash
 $VSEARCH --cluster_size derep.fasta \
     --threads $THREADS \
     --id 0.97 \
     --strand plus \
     --sizein \
     --sizeout \
     --fasta_width 0 \
     --centroids centroids.fasta
```
---
### Sort centroids and remove singletons.
This sorts the sequences by abundance and removes the sequences found only once
 ```bash
 $VSEARCH --sortbysize centroids.fasta \
     --threads $THREADS \
     --sizein \
     --sizeout \
     --fasta_width 0 \
     --minsize 2 \
     --output sorted.fasta
```
---
### denovo chimera detection.
Looks for and removes chimeras formed from clustering
``` bash
$VSEARCH --uchime_denovo sorted.fasta \
     --sizein \
     --sizeout \
     --fasta_width 0 \
     --qmask none \
     --nonchimeras denovo.nonchimeras.fasta
```
---
### Reference chimera detection against SILVA138.
Blasts sequences against SILVA138 database. Removes chimeras
 ``` bash
 $VSEARCH --uchime_ref denovo.nonchimeras.fasta \
     --threads $THREADS \
     --db silvamod138pr2.fasta \
     --sizein \
     --sizeout \
     --fasta_width 0 \
     --qmask none \
     --dbmask none \
     --nonchimeras nonchimeras.fasta
```
---
### Relabel OTUs.
 ``` bash
 $VSEARCH --fastx_filter nonchimeras.fasta \
     --threads $THREADS \
     --sizein \
     --fasta_width 0 \
     --relabel OTU_ \
     --fastaout XXX_otus.fasta
```
 ---
### map sequences to OTU
Change "XX" to unique file name for dataset
 ``` bash
 $VSEARCH --usearch_global all.fasta \
     --threads $THREADS \
     --db XXX_otus.fasta \
     --id 0.97 \
     --strand plus \
     --sizein \
     --sizeout \
     --fasta_width 0 \
     --qmask none \
     --dbmask none \
     --otutabout XXX.otutab.txt 
```
---
#### Sort OTU table numerically 
 ``` bash
 sort -k1.5n XXX.otutab.txt > XXX.otutab.sorted.txt
 ```
---

## LULU cleanup
make blast db in terminal using last .fasta file from "relabel OTUs" step
``` bash
makeblastdb -in Dec22exp_otus.fasta -dbtype nucl

blastn -db Dec22exp_otus.fasta -outfmt '6 qseqid sseqid pident' -out Dec22exp_LULU_match_list.txt -qcov_hsp_perc 95 -perc_identity 95 -query Dec22exp_otus.fasta -num_threads 3
```

Open R Studio

Run LULU.R script and change "XX" to unique name used in output from mapping sequences to OTUs 
``` bash
setwd("/home/hba052/Documents/University/Sequencing/data/Deep_biosphere/XX")

require(lulu)
require(methods)

matchlist = read.delim("XX_LULU_match_list.txt", header=FALSE,as.is=TRUE, stringsAsFactors=FALSE)
otus.all = read.delim("XX_all.otutab.sorted.txt",row.names=1,header=T,sep="\t")
curated_result <- lulu(otus.all,matchlist, minimum_match = 97)
lulus = curated_result$curated_table
write.table(data.frame("OTU"=rownames(lulus),lulus),"XX_table_curated.tsv", row.names=FALSE, 
            quote=F, sep="\t")

  curated_result$curated_table 

```

### Remove sequences from fasta file that got removed using LULU
(in Python)
``` bash
python filter_seq_by_OTUtable.py Dec22exp_otus.fasta Dec22exp_table_curated.tsv > Dec22exp_OTUs_curated.fasta
```
________________________________________________________________________

## CREST Classifying
Using CREST4 which is installed locally on my computer:
activateconda environment to run CREST

``` bash
conda activate condaHRB![[Orcutt_VSEARCH_outputs.xlsx]]
```

Run CREST4
``` bash
crest4 -f Dec22exp_OTUs_curated.fasta
```

Rename the OTU files to have unique name for datasets you are working on
"**XX_OTU_table.csv "

-----------------------------------------------------------------------
												 IF USING KRAFTLYNG 

Make xml file
``` bash
blastn -query OTUs_curated.fasta -db /home/hba052/Documents/University/Sequencing/####working file####/silvamod138pr2.fasta -max_target_seqs 100 -outfmt 5 -out otus.xml -num_threads 6
```
you might need to upload the silvamod138 data base to kraftlyng
_______________________________________________________________________

## Merge assignments and OTUs
In excel or spreadsheets, merge the OTU table with the assignments from the CREST classifying. Save as your OTU table is .csv format

Carefully look through this for errors. 