# PolSter
Pol II density estimated by statistical inference of transcription elongation rate by total RNA-seq
## Overview
The RNA-seq library was prepared from rRNA-depleted non-poly-A transcripts (total RNA-seq) provides a transcriptomic profile of nascent RNAs undergoing transcription with co-transcriptional splicing. In general, the RNA-seq reads exhibit a sawtooth pattern in a gene, which is characterized by a monotonically decreasing gradient across introns in the 5’ to 3’ direction, and by substantially higher levels of RNA-seq reads present in exonic regions. Such patterns result from the process of underlying transcription elongation by RNA polymerase II (Pol II). The objective is to reconstruct the spatial distribution of transcription elongation rates in a gene from a given noisy, sawtooth-like profile.


## Demonstoration
Run example Total RNA-seq dataset [1] (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE36799) of this project.
These were the following commands that were used to estimate Pol II density in mouse dataset.

1. Download input data "input_data.zip"
2. gcc -lm Estimate_Pol2.c
3. ./a.out /DIRECTORY_input_data/

To interpret transcription elongation rate, you can create such figures as "Estimated_result_example.jpg". The Pol II existence probability is inversely proportional to the elongation rate.

## Usage

### Input formats
Accepted input formats for this software include txt files. For example, with a txt input, in order to build models, we expect the following columns to have entries for each row in the txt file. Each intron was divided into bins with intervals equal to 400 bp. An exonic region was treated positionally as a single point.
Here's an example of a potential input file.

*_read.txt: read count data across genomic bin positions from TSS, * shows gene number. 

*_EI.txt: exon or intron across genomic bin positions, * shows gene number

Particle_num.txt: The number of particles for particle filter: recommend more than 100000

SRR960177_pickup_gene_num_joint_ei_av_mu_tau_sigma.txt: hyper parameter for each gene, mu is the initial state of variable of state model, tau is noise of system model, and sigma is noise of meaturement model (log). The detail is described in our paper as follows. 

SRR960177_pickup_gene_num_joint_ei_av_chrPosSE_sig_len5cor05_cov01raw.txt: gene number, chrosome numer, gene name, strand, start, end. Here, one gene resulted in combined all isoforms in one.

## Reference

[1] Sigova AA, Mullen AC, Molinie B, Gupta S et al., Divergent transcription of long noncoding RNA/mRNA gene pairs in embryonic stem cells., *Proc Natl Acad Sci U S A*, 2013 Feb 19;110(8):2876-81., PMID: 23382218


## Credit

If you use this program in your work, please cite:

Yumi Kawamura, Shinsuke Koyama, and Ryo Yoshida., Statistical inference of the rate of RNA polymerase II elongation by total RNA sequencing., *Bioinformatics*, 2018, *in preparation*.

