# HAPpy-ABCENTH
A pipeline for homology based annotation utilizing exon structural conservation

## Installation/dependencies
Dependencies: python (inlcuding numpy, ete3, and intervaltree libraries), HMMER version 3, MAFFT. Wise2 is also required to use the "full_pseudoexon_search" mode in ABCENTH. All of these dependencies can be installed using conda (see below).
### Install with conda + pip
```
conda create -n happy -c bioconda -c etetoolkit python=3 numpy ete3 wise2 mafft=7 hmmer=3 intervaltree
conda activate happy
pip install HAPpy-ABCENTH
```
### Install from source
```
git clone https://github.com/biorover/HAPpy-ABCENTH.git
cd HAPpy-ABCENTH
conda env create -n happy -f conda_env.yaml
conda activate happy
python setup.py install
```
## Usage
HAPpy can be used in many different ways (full options list available with `HAPpy -h`), but there's two particularly common use cases.
1. Generic homology-based gene structure prediction (note: requires genewise, also many other tools currently outperform HAPpy for this):
    `HAPpy --genome <target genome fasta> --protein_seqs <protein fasta> [--threads <nthreads>]`
 
 2. Annotation of tandemly duplicated genes from large gene families with conserved exon structure (e.g. insect odorant receptors)
    If you have annotations in gff form:
    ```
    HAPpy --genome <target genome fasta> --ref_genome <one or more reference genome fastas> \
    --annotations <one gtf per ref genome> --cutoff <p distance at which proteins are clustered, 0.45 for ant ORs> \
    --search_mode exons --annotator ABCENTH
    ```
    If you have pre-built hmm files for each exon in each cluster:
    ```
    HAPpy --genome <target genome fasta> --hmm_dir <path to hmm folder> \
    --annotator ABCENTH --threads <# threads> --output_dir <output directory>
    ```
    Either way, the output gtf will be in \<output directory\>/ABCENTH.gtf
    If a gene name ends in a number, it's a complete model. Otherwise, the following letters indicate gene model issues:
        - P: Pseudogene (no compatible splice sites for exon, frameshift in exon, or stop codon in exon)
        - N: Missing N terminal
        - C: Missing C terminal
        - I: Missing internal section

ABCENTH can be used directly given a correctly formated table of exons hits and table of exon structures for each gene (will write up file specifications soon).
Run `ABCENTH -h` for full list of options.

## The pipeline
This package consists of two separate programs: HAPpy (Homology Annotation Pipeline) and ABCENTH (Annotation Based on Conserved Exons Noticed Through Homology). HAPpy locates putative genetic loci through homology and feeds them into a gene predictor, which can be ABCENTH (also can use genewise).

### HAPpy
The HAPpy pipeline is tailored to finding the loci for highly divergent gene expansions. In it's simplest form, HAP will take input sequences (either as HMMER hmm files, protein fasta files, or gtfs with a reference genome), build HMM profiles (if the input isn't already an HMM profile), and match these against all orfs in a target genome using hmmsearch (via the wrapper thammerin). It can then parse the target genome into candidate loci and feed these into genewise. Overall, not very exciting- but theres two optional settings that do more interesting things.

  --cutoff <float 0.0-1.0> : if a value less than 1 (and greater than 0) is given, the reference protein sequences are clustered based on pairwise similarity (using UPGMA clustering). This is mostly useful for large multi-gene families with high divergence, as HAP will then build HMMs for each cluster which will significantly boost detection sensitivity (and save time since fewer total HMMs will need to be aligned to the target genome).
  
  --search_mode <"fl" or "exons"> : The default "fl" (for "full length") mode uses the entire protein sequence to find candidate exons, but if you provide reference sequences in the form of reference genomes and gtf files (--ref_genome and --annotations) and set --search_mode to "exons", HAP.py will build HMMs seperately for each exon. When exon structures are conserved, this can increase sensitivity and precision. Furthermore, it allows you to use the ABCENTH annotator instead of Genewise (if you expect exon phase to be conserved).
  
For some gene families (e.g. ant odorant receptors, for which this program was originally written), exon structure is perfectly conserved for all genes that cluster at a specific cutoff (about 0.5-0.45 for ant ORs), thus using both "--cutoff <float>" and "--search_mode exons" is an extremely powerful way to find all exons in a genome.
   
### ABCENTH
ABCENTH is a gene finder specifically devised for multigene families with extremely high sequence divergence but highly conserved exon structure. It is also designed to avoid gene fusion in tandem arrays. It works by taking candidate exons from a table of homology search hits and extending the candidate region by finding splice sites that recapitulate the expected exon length and phase. Phase is completely constrained, while length can be variable- the program will first look for splice sites at the exact right positions and then will alternately look at +3,-3,+6,-6, etc (until +30,-30, although this limit can be edited by changing the default "max_offset" in the exon_finder function at line 57) until it finds a valid splice site. If no valid splice site is found, the exon is still returned but is flagged as pseudogenized. As each exon is discovered, they are strung together into genes assuming they fall into the correct order (e.g. exon 2 follows exon 1). If an exon is skipped (e.g. a candidate exon 3 is found after a candidate exon 1, or a candidate exon 2 is found with no candidate exon 1), ABCENTH will extract all possible exons (all DNA segments between AG and GT which have an ORF with the correct phase and are roughly the right lenght) and select the one with the best HMMER match to the missing exon. It is thus highly robust to short exons that are often missed by homology searches (provided of course the phase and rough length of the exon is conserved). To run, it therefore needs a table of candidate hits, information about the expected length, phase, and position within the gene of each exon query, and HMM profiles for each exon query (against which to test candidate missed exons). These are all automatically generated by HAPpy, so I won't document them here unless someone requests it under "issues" (assumedly because they want to run ABCENTH without HAPpy).
