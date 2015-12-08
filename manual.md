**USAGE**:
orthAgogue -i 

<file\_path>

 [-t][-p][-s][-O][-P][-S][-A][-w][-b][-c][-u|-e][-o][-C][m](m.md)
<br>

<b>INPUT</b>:<br>
<br>
-i    --input                Path to the BLAST output file (string, mandatory). The file must be in the tabular '-m 8' format (12 fields, one row per HSP (high scoring pair).<br>
<br>
-t    --taxon_index          Index of the field containing taxon IDs in columns 1 and 2 (positive integer, default '0').<br>
<br>
-p    --protein_index        Index of the field containing protein IDs in columns 1 and 2  (positive integer, default  '1').<br>
<br>
-s    --seperator            Field separator used in columns 1 and 2 (only non-integers and -non-decimals, default the pipe '|').<br>
<br>
<b>OUTPUT</b>:<br>
<br>
-O    --output_dir           Output directory path (string, default current directory).<br>
<br>
-P    --pipe                 Send the output in the form of MCL native matrix (all.mcl) to STDOUT for piping.<br>
<br>
-S    --sort                 Sort .abc files by the scores.<br>
<br>
<br>
<b>OPERATIONAL</b>:<br>
<br>
-c    --cpu                  The numbers of threads to use (integer, default '2').<br>
<br>
-u    --use_scores           Use BLAST scores (column 12) instead of e-values (default, column 11), mutually exclusive with -e.<br>
<br>
-A    --all_to_norm          Use all protein pairs retained after filtering to compute normalization basis (by default only protein pairs forming homology relations).<br>
<br>
-w    --without_norm         Skip normalization of similarity scores.<br>
<br>
-b    --best_hsp             In case of multiple HSPs for a protein pair use only the best one (for emulating OrthoMCL only).<br>
<br>
<br>
<b>FILTERING</b>:<br>
<br>
-e    --cutoff               The e-value cutoff (positive floating number, see examples below), mutually exclusive with -u.<br>
<br>
-o    --overlap              Threshold for protein pair overlap (integer [1-100], see examples below).<br>
<br>
-C    --strict_coorthologs   Restricted definition of co-orthologs, relations between inparalogs of orthologs are excluded. Useful to give more weight to orthologs during MCL computation.<br>
<br>
-m    --min_proteins         Threshold for the number of proteins in a proteome. Useful for handling files containing many taxa with just a few proteins, e.g. the complete SwissProt.<br>
<br>

<b>EXAMPLES</b>
<br>
1. orthAgogue -i myblast.out -e 6 -O myoutdir # excludes protein pairs with e-values above 1e-06 and saves results in 'myoutdir'.<br>
<br>
2. orthAgogue -i myblast.out -e 6 -o 50 -O myoutdir # the same as above but excluding protein pairs with the overlap below 50%.<br>
<br>
3. orthAgogue -i myblast.out -u -o 50 -O myoutdir # the same as above but without filtering by e-values and using BLAST scores instead of e-values in order to resolve HSPs with the '0.0' e-value; the required e-value cutoff should be set while running BLAST.<br>
<br>
4. orthAgogue -i myblast.out -b -e 6 -O myoutdir # the same as 1 but using only the best HSP for any pair of proteins (OrthoMCL emulation).<br>
<br>

See also <a href='https://code.google.com/p/orthagogue/'>https://code.google.com/p/orthagogue/</a>
<br><br>
The software was developed by O.K. Ekseth under supervison of Dr. V.N. Mironov. Questions to be forwarded to orthagogue-issue-tracker@googlegroups.com