# Patera
Scripts for a study on the status of Patera nantahala

noGaps_nucleotides.sh is a rather complicated AWK command that removed alignment positions with all Ns. This is useful for RADseq data when some loci/contigs do not have sequence overlaps. As written, fasta files must have <.fas> extension and sequences must be on a single line. 

niche-modeling_Patera_FINAL.R is the script used for environmental niche modeling. See comments in script, but this will take a long time to run, needs lots of memory, and you may want to split commands across multiple R scripts to speed things up.

Patera_SNAQ.jl script is in Julia. It is the code I used to Run SNAQ.
Patera_SNAQ_plotting.R plots some output from SNAQ
