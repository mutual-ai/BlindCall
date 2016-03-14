# BlindCall

ultra-fast base-calling of high-throughput sequencing data by blind deconvolution.

Ye C, Hsiao C, Corrada Bravo H

http://bioinformatics.oxfordjournals.org/content/30/9/1214.long


Run_201305_PhiX174_100nt_pe_UMD_HiSeq.m is the main file. 

Test data can be found:
ftp://ftp.cbcb.umd.edu/pub/data/hcorrada/BlindCall_data.tar.gz

To run specific datasets, 
1) Modify the directory 'CIF_dir' to the directory containing the cif files.
2) Modify the 'totCycles' to the total cycles of a run.
3) If the cif directory contains data for several tiles, modify 'Tiles_beg' and 'Tiles_end' to specify which tiles to run.
