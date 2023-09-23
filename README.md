
# Bulk_metagenome_virome_2023

Bioinformatics workflow used in the article:

> Sabatino R, Sbaffi T, Corno G, Fontaneto D, Periyasamy S, Di Cesare A. 2022.  Bacteriophages limitedly contribute to the antimicrobial resistome of microbial communities in wastewater treatment plants. Microbiology Spectrum,e01101-23. doi: [xxxx](https://doi.org/10.1128/spectrum.01101-23
PDF/EPUB

).

For the biostatistics workflow of the same paper, [Click here](Intra_Extra_DNA_script_statistical_analysis.R)


## Contacts

**Tomasa Sbaffi**  
Postdoctoral Researcher, bioinformatics  
[E-mail](mailto:tomasa.sbaffi@gmail.com)

**Andrea Di Cesare**  
Principal Investigator  
[E-mail](mailto:andrea.dicesare@cnr.it)


## Table of contents

1. [Before starting](#before-starting)
2. [Read based analysis](#read-based-analysis)
3. [Assembly based analysis](#assembly-based-analysis)


## Before starting

### You will need to have these softwares installed and added to your path

* SRA Toolkit v2.11.3: https://github.com/ncbi/sra-tools
* fasterq-dump v2.10.8: https://github.com/glarue/fasterq_dump
* GNU parallel: https://www.gnu.org/software/parallel
* BLAST v2.10.1: https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+
* seqtk v1.3: https://github.com/lh3/seqtk
* CoverM v0.6.1: https://github.com/wwood/CoverM/
* fastQC v0.11.9: http://www.bioinformatics.babraham.ac.uk/projects/fastqc/
* multiQC v1.8: https://multiqc.info/
* Cutadapt v1.16: https://cutadapt.readthedocs.io/
* MEGAHIT v1.1.2.9: https://github.com/voutcn/megahit/
* bowtie v2.3.5: http://bowtie-bio.sourceforge.net/bowtie2/
* SAMtools v1.9: http://www.htslib.org/
* GTDB-Tk v2.1.0 and GTDB release 05-RS95: https://gtdb.ecogenomic.org/
* MetaPhlan v3.0.14
* Deeparg v1.0.2
* Trimmomatic v0.39
* vsearch v2.17.1
* Trimgalore v0.6.5
* bwa-mem v0.7.17
* MetaBAT v2.15.0
* RefineM v0.1.2
* CheckM v1.2.0
* Magpurify v2.1.2
* prodigal v2.6.3
 
It should be possible running the analysis with any UNIX-based OS.

### Define number of threads to use

You should change below to the number of cores available in your system and set you project name:

```bash
NTHREADS=40
project=THIS_PROJECT
```


### Note on the data used for this workflow

Eight coassemblies of two samples were used, this work is a follow up of: 

>Periyasamy S, Sabbatino R, Sbaffi T, Fontaneto D, Corno G, Di Casare A. 2022. Extracellular DNA includes an important fraction of high-risk antibiotic resistance genes in treated wastewaters. Environmental Pollution

The data here used is assembled as described in the [WWTP Intra_extra_cellular DNA metagenomics 2022](https://github.com/TomasaSbaffi/WWTP_extracell_DNA_metagenomics_2022#wwtp-intra_extra_cellular-dna-metagenomics-2022) project.


## 1. Binning, with vamb

Singularly for each coassembly, using the input from megahit and the depth files, vamb calculated bins of contigs, the minimum size of putative bins was set to 3500.

```bash
for object in IVB_pre IVB_post EVB_pre EVB_post ICB_pre ICB_post ECB_pre ECB_post
do
vamb --outdir ${object}/vamb/ --fasta ${object}_metasens_m1000/final.contigs.fa --minfasta 3500 --jgi ${object}/depth.txt # minfasta set to 3500
done
```


## 2. PHAMB pipeline

PHAMB pipeline to discerne viral bins from non viral bins. This process was run singularly for each coassembly.
The PHAMB pipeline consists in the following steps:

- DVF (DeepVirFinder) was run on final.contigs.fa (-l 1000 -c 1) using its own database which (check this!) comes from UNIProt-Virus

```bash
conda activate dvf_env

dvf_path="/path/to/DeepVirFinder/database"

ln -s ../../IVB_pre_metasens_m1000/final.contigs.fa .  #correct with loop here

mkdir annotations

python3 ${dvf_path}/dvf.py -i final.contigs.fa -o DVF -l 1000 -c 1
mv DVF/final.contigs.fna_gt2000bp_dvfpred.txt annotations/all.DVF.predictions.txt
```

- prodigal (-p meta -g 11) was used to predict ORFs from final.contigs.fa
- hmmsearch (-E 1.0e-5) was used to align coassemblies ORFs to the MicompleteDB
- hmmsearch (-E 1.0e-5) was used to align coassemblies ORFs to the VOGDB

```bash
conda activate phamb_env

micompleteDB="/path/to/Bact105.hmm"
VOGDB="/path/to/AllVOG.hmm"

for object in IVB_pre IVB_post EVB_pre EVB_post ICB_pre ICB_post ECB_pre ECB_post
do
phamb_path="/path/to/${object}/PHAMB_out"

prodigal -i ${phamb_path}/final.contigs.fa -d ${phamb_path}/genes.fna -a ${phamb_path}/proteins.faa -p meta -g 11

hmmsearch --cpu 40 -E 1.0e-05 -o ${phamb_path}/output1.txt --tblout ${phamb_path}/annotations/all.hmmMiComplete105.tbl ${micompleteDB} ${phamb_path}/proteins.faa

hmmsearch --cpu 40 -E 1.0e-05 -o ${phamb_path}/output2.txt --tblout ${phamb_path}/annotations/all.hmmVOG.tbl ${VOGDB} ${phamb_path}/proteins.faa
done
```

