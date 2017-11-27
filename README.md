# Removing-mouse-contamination-from-WES-
The major issue of The WES data from patient derived xenofraft samples (i.e, human tissue that is engrafted in mouse) is mouse contamination. This program removes mouse reads and ambiguous reads from the aligned BAM files. This program was written for the BAM  files aligned with BWA-MEM.

The program take three  arguments:  mouseAlignment humanAlignment outputPrefix

The program has 10 outputs : Since this is a paired end reads : 

Human_only reads  = Aligned perfectly to human reference
Human_better reads = Both the pair has better score for human alignment 
Human_better_one reads = one of the pair has better score for human alignment 
Human_better_mixedreads = Difference in  pair score favouring  human alignment
Mouse_only reads= Aligned perfectly to mouse reference
Mouse_better reads= Both the pair has better score for mouse alignment 
Mouse_better_one reads= one of the pair has better score for mouse alignment 
Mouse_better_mixed reads=Difference in pair score favouring mouse alignment
Mixed reads = Scores are mixed in pairs cannot decide
Tie reads = scores are tied.


