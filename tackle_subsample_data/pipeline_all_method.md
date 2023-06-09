# How we run all pipeline in our paper
## Basecall and alignment
```sh
## basecall and alignment was applied on all dataset
guppy_basecaller -i ./single/ -s ./guppy_out -c rna_r9.4.1_70bps_hac.cfg --device auto --recursive
cat */*.fastq>all.fastq
minimap2 -ax map-ont -t --MD16 SINV_Toto1101.fa all.fastq | samtools view -hbS -F 260 - | samtools sort -@ 6 -o all.bam
```
## Quality control
we used the tool from ONT called `pomoxis`, and plot the figure by our own code
```sh
##
stats_from_bam all.bam >stats_filtered.txt
```

## Epinano
EpiNano_Error utilizes the differences between two alignment files in basecalling errors (mismatches, insertions, and deletions) and the alterations in per-base qualities to predict modified bases.
###  mismatch
```sh
python path/EpiNano/Epinano_Variants.py -n 20 -R SINV_Toto1101.fa -b ivt.bam -s path/EpiNano/misc/sam2tsv.jam --type t
python path/EpiNano/Epinano_Variants.py -n 20 -R SINV_Toto1101.fa -b wt.bam -s path/EpiNano/misc/sam2tsv.jam --type t
Rscript path/EpiNano/Epinano_DiffErr.R -k ivt/ivt.plus_strand.per.site.csv -w wt/wt.plus_strand.per.site.csv -o epinano_mismatch_ -f mis -d 0.1 -p
```
###  sumErr
```sh
python path/EpiNano/misc/Epinano_sumErr.py --file wt.plus_strand.per.site.csv --out wt.sum_err.csv --kmer 0
python path/EpiNano/misc/Epinano_sumErr.py --file ivt.plus_strand.per.site.csv --out ivt.sum_err.csv --kmer 0
Rscript path/EpiNano/Epinano_DiffErr.R -k ivt/ivt.sum_err.csv -w wt/wt.sum_err.csv -o epinano_sumErr -f sum_err -d 0.1 -p
```
## Differr
```sh
python main.py -b wt/wt.bam -a ivt.bam -r SINV_Toto1101.fa -o differr.bed
```
## ELIGOS2
```sh
 eligos2 pair_diff_mod -tbam wt.bam  -cbam ivt.bam  -reg gene.bed -ref SINV_Toto1101.fa -t 16 --pval 1 --oddR 0 --esb 0 -o eligos2_results 
```
## Nanocompore
```sh
nanopolish index -d single/ all.fastq
nanopolish eventalign --reads all.fastq --bam all.bam --genome SINV_Toto1101.fa --print-read-names --scale-events --samples > eventalign_reads.tsv
nanocompore eventalign_collapse -t 6 -i eventalign_reads.tsv -o eventalign_collapsed_reads_tsv
nanocompore sampcomp \
    --file_list1 ./wt/eventalign_collapsed_reads_tsv/out_eventalign_collapse.tsv 
    --file_list2 ./ivt/eventalign_collapsed_reads_tsv/out_eventalign_collapse.tsv 
    --label1 wt
    --label2 ivt
    --fasta SINV_Toto1101.fa \
    --outpath ./nanocompore_results
```
## Tombo
### Tombo_comp
```sh
# Single format is required to run tombo, tun ont-fast5-api firstly to transfer if needed.
tombo preprocess annotate_raw_with_fastqs --fast5-basedir single/ --fastq-filenames all.fastq
tombo resquiggle single/ SINV_Toto1101.fa --processes 20 --num-most-common-errors 5 --overwrite
# level_sample_compare 
tombo detect_modifications level_sample_compare --fast5-basedirs wt/single/ \
    --alternate-fast5-basedirs ivt/single/ \
    --statistics-file-basename sample.level_compare_sample
    
tombo text_output browser_files --statistics-filename sample.level_compare_sample.stats\
    --browser-file-basename sample.level_samp_comp_detect --file-types statistic
    
tombo text_output browser_files --statistics-filename sample.level_samp_comp_detect.tombo.stats \
  --fast5-basedirs wt/single/ --control-fast5-basedirs ivt/single/ \
  --browser-file-basename sample.level_samp_comp_detect --file-types difference
```
### Tombo_de novo
```sh
tombo detect_modifications de_novo --fast5-basedirs single/ \
    --statistics-file-basename sample.de_novo
```
## Xpore
```sh
nanopolish eventalign  --reads ivt.fastq --bam ivt.bam  --genome SINV_Toto1101.fa --scale-events > ivt.eventalign.txt
nanopolish eventalign  --reads wt.fastq --bam wt.bam  --genome SINV_Toto1101.fa --scale-events > wt.eventalign.txt

xpore dataprep --eventalign ivt.eventalign.txt --out_dir ivt_dataprep
xpore dataprep --eventalign wt.eventalign.txt --out_dir wt_dataprep

### create a yml file
data:
    IVT:
        rep1: ./IVT/IVT_dataprep/
    wt:
        rep1: ./wt/wt_dataprep/

out: ./xpore_results # output dir
xpore diffmod --config xpore.yml
```

## DRUMMER
```sh
python DRUMMER.py -r SINV_Toto1101.fa -t ivt.bam -c wt.bam -p 1 -o DURMMER_result -a exome -n NC_001547.1 -m True
```
## m6anet
```sh
nanopolish eventalign --reads all.fastq --bam all.bam --genome SINV_Toto1101.fa \
--scale-events --signal-index --summary /path/to/summary.txt  --threads 50 > eventalign.txt

m6anet dataprep --eventalign eventalign.txt \
                --out_dir m6anet_dataprep  --n_processes 4
                
m6anet inference --input_dir m6anet_dataprep --out_dir m6anet_output  --n_processes 4 --num_iterations 1000
```

## Nanom6A
```sh
multi_to_single_fast5 -i guppy -s single -t 40 --recursive

tombo resquiggle --overwrite --basecall-group Basecall_1D_001 single/ SINV_Toto1101.fa  \
--processes 40 --fit-global-scale --include-event-stdev

find single -name "*.fast5" >files.txt

extract_raw_and_feature_fast --cpu=20 --fl=files.txt -o result --clip=10

predict_sites --cpu 20 -i result -o result_final -r SINV_Toto1101.fa -g SINV_Toto1101.fa -b gene2transcripts.txt 
```