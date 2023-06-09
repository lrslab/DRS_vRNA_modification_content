# How we run all pipeline in our paper
## Basecall and alignment
```sh
guppy_basecaller -i ./single/ -s ./guppy_out -c rna_r9.4.1_70bps_hac.cfg --device auto --recursive
cat */*.fastq>final.fastq
minimap2 -ax map-ont -t --MD16 SINV_Toto1101.fa final.fastq | samtools view -hbS -F 260 - | samtools sort -@ 6 -o xxx.bam
```
## Quality control
we used the tool from ONT called `pomoxis`
```sh
##
stats_from_bam xxx.bam >stats_filtered.txt
```

## Epinano
EpiNano_Error utilizes the differences between two alignment files in basecalling errors (mismatches, insertions, and deletions) and the alterations in per-base qualities to predict modified bases.
###  Mismatch
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
## differr
```sh
python main.py -b wt/wt.bam -a ivt.bam -r SINV_Toto1101.fa -o differr.bed
```
## nanocompore
```sh
nanopolish index -s {sequencing_summary.txt} -d {raw_fast5_dir} {basecalled_fastq}
nanopolish eventalign --reads {basecalled_fastq} --bam {aligned_reads_bam} --genome {transcriptome_fasta} --print-read-names --scale-events --samples > {eventalign_reads_tsv}
nanocompore eventalign_collapse -t 6 -i {eventalign_reads_tsv} -o {eventalign_collapsed_reads_tsv}
nanocompore sampcomp \
    --file_list1 ./data/S1_R1.tsv
    --file_list2 ./data/S2_R1.tsv
    --label1 ivt
    --label2 wt
    --fasta ./reference/ref.fa \
    --outpath ./results
```
## tombo
```sh
# Single format is required to run tombo, tun ont-fast5-api firstly to transfer if needed.
tombo preprocess annotate_raw_with_fastqs --fast5-basedir <fast5s-base-directory> --fastq-filenames all.fastq
tombo resquiggle <fast5s-base-directory> SINV_Toto1101.fa --processes 20 --num-most-common-errors 5 --overwrite
# level_sample_compare 
tombo detect_modifications level_sample_compare --fast5-basedirs <fast5s-base-directory> \
    --alternate-fast5-basedirs <alternate-fast5s-base-directory> \
    --statistics-file-basename sample.level_compare_sample
tombo text_output browser_files --statistics-filename sample.level_compare_sample.stats\
    --browser-file-basename sample.level_samp_comp_detect --file-types statistic
tombo text_output browser_files --statistics-filename sample.level_samp_comp_detect.tombo.stats
  --fast5-basedirs <fast5s-base-directory> --control-fast5-basedirs <control-fast5s-base-directory> --browser-file-basename sample.level_samp_comp_detect --file-types difference
```
## xpore
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

