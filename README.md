# Anopheles-assembly
A project conducted during the work in the institute of cytology and genetics of SB RAS, laboratory of delelopmental genetics, group leader — Veniamin Fishman, PhD.

Existed genome assembleys of mosquitoes from the Anopheles genus were incomplete. They had large empty spaces, which couldn't be filled with the use of the second generation sequencing.
Here, we used the data of the third generation sequencing (Nanopore, PacBio) to fill these gaps.

The script uses the pysam module to work with the *.sam file with long reads. It searches for reads overlapping empty regions of the genome. Once found, the read is stored with its metadata in *.fa format. Futher, multiple sequence alignment
has to be performed with these reads to define the consensus sequence and then it is integrated into the genome where the empty space has been.
