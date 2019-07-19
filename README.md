# HAPpy-ABCENTH
A pipeline for homology based annotation utilizing exon structural conservation

## Installation/dependencies
This is a ready-to-go python script- as long as all of the python modules are left in the same directory it should work (provided the dependencies are installed).

Dependencies: python (inlcuding numpy, ete2, and a few other common libraries), HMMER version 3, MAFFT, thmmerin.py HMMER wrapper (included in git archive but needs to be in the system path)

## Usage
HAPpy can be used in many different ways (full options list available with "python HAP.py -h"), but there's two particularly common use cases.
1. Generic homology-based gene structure prediction (note: requires genewise, also many other tools currently outperform HAPpy for this):
    python HAP.py --genome <target genome fasta> --protein_seqs <protein fasta> [--threads <nthreads>]
 
 2. Annotation of tandemly duplicated genes from large gene families with conserved exon structure (e.g. insect odorant receptors)
    python HAP.py --genome <target genome fasta> --ref_genome <one or more reference genome fastas> --annotations <one gtf per ref genome> --cutoff <p distance at which proteins are clustered, 0.5 for ant ORs> --search_mode exons --annotator ABCENTH

## The pipeline
This package consists of two separate programs: HAP.py (Homology Annotation Pipeline) and ABCENTH (Annotation Based on Conserved Exons Noticed Through Homology). HAP.py locates putative genetic loci through homology and feeds them into a gene predictor, which can be ABCENTH (also can use genewise).

### HAP.py
The HAP.py pipeline is tailored to finding the loci for highly divergent. In it's simplest form, HAP will take input sequences (either as HMMER hmm files, protein fasta files, or gtfs with a reference genome), build HMM profiles (if the input isn't already an HMM profile), and match these against all orfs in a target genome using hmmsearch (via the wrapper thmmerin.py). It can then parse the target genome into candidate loci and feed these into genewise. Overall, not very exciting- but theres two optional settings that do more interesting things.

  --cutoff <float 0.0-1.0> : if a less than 1 (and greater than 0) is given, the reference protein sequences are clustered based on pairwise similarity (using UPGMA clustering). This is mostly useful for large multi-gene families with high divergence, as HAP will then build HMMs for each cluster which will significantly boost detection sensitivity (and save time since fewer total HMMs will need to be aligned to the target genome).
  
  --search_mode <"fl" or "exons"> : The default "fl" (for "full length") mode uses the entire protein sequence to find candidate exons, but if you provide reference sequences in the form of reference genomes and gtf files (--ref_genome and --annotations) and set --search_mode to "exons", HAP.py will build HMMs seperately for each exon. When exon structures are conserved, this can increase sensitivity and precision. Furthermore, it allows you to use the ABCENTH annotator instead of Genewise (if you expect exon phase to be conserved).
  
For some gene families (e.g. ant odorant receptors, for which this program was originally written), exon structure is perfectly conserved for all genes that cluster at a specific cutoff (about 0.5-0.45 for ant ORs), thus using both "--cutoff <float>" and "--search_mode exons" is an extremely powerful way to find all exons in a genome.
   
### ABCENTH
Explanation coming
