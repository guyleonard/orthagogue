# orthAgogue's Program Interface #
We include here how orthAgogue is used; simply launch it without any arguments. (In practical terms: _orthAgogue_ or _./orthAgogue_.) In the below encapsulated description, we've included a synopsis of usage, a list of options explained and a number of examples of use.

```bash


USAGE:
orthAgogue -i <file_path> [-t][-p][-s][-O][-P][-S][-A][-w][-b][-c][-u|-e][-o][-C][m]

INPUT:
-i    --input                Path to the BLAST output file (string, mandatory). The file must be in the tabular '-m 8' format (12 fields, one row per HSP (high scoring pair).
-t    --taxon_index          Index of the field containing taxon IDs in columns 1 and 2 (positive integer, default '0').
-p    --protein_index        Index of the field containing protein IDs in columns 1 and 2  (positive integer, default  '1').
-s    --seperator            Field separator used in columns 1 and 2 (only non-integers and -non-decimals, default the pipe '|').

OUTPUT:
-O    --output_dir           Output directory path (string, default current directory).
-P    --pipe                 Send the output in the form of MCL native matrix (all.mcl) to STDOUT for piping.
-S    --sort                 Sort *.abc files by the scores.

OPERATIONAL:
-c    --cpu                  The numbers of threads to use (integer, default '2').
-u    --use_scores           Use BLAST scores (column 12) instead of e-values (default, column 11), mutually exclusive with -e.
-A    --all_to_norm          Use all protein pairs retained after filtering to compute normalization basis (by default only protein pairs forming homology relations).
-w    --without_norm         Skip normalization of similarity scores.
-b    --best_hsp             In case of multiple HSPs for a protein pair use only the best one (for emulating OrthoMCL only).

FILTERING:
-e    --cutoff               The e-value cutoff (positive floating number, see examples below), mutually exclusive with -u.
-o    --overlap              Threshold for protein pair overlap (integer [1-100], see examples below).
-C    --strict_coorthologs   Restricted definition of co-orthologs, relations between inparalogs of orthologs are excluded. Useful to give more weight to orthologs during MCL computation.
-m    --min_proteins         Threshold for the number of proteins in a proteome. Useful for handling files containing many taxa with just a few proteins, e.g. the complete SwissProt.

EXAMPLES
1. orthAgogue -i myblast.out -e 6 -O myoutdir # excludes protein pairs with e-values above 1e-06 and saves results in 'myoutdir'.
2. orthAgogue -i myblast.out -e 6 -o 50 -O myoutdir # the same as above but excluding protein pairs with the overlap below 50%.
3. orthAgogue -i myblast.out -u -o 50 -O myoutdir # the same as above but without filtering by e-values and using BLAST scores instead of e-values in order to resolve HSPs with the '0.0' e-value; the required e-value cutoff should be set while running BLAST.
4. orthAgogue -i myblast.out -b -e 6 -O myoutdir # the same as 1 but using only the best HSP for any pair of proteins (OrthoMCL emulation).

See also https://code.google.com/p/orthagogue/

The software was developed by O.K. Ekseth under supervison of Dr. V.N. Mironov. Questions to be forwarded to orthagogue-issue-tracker@googlegroups.com
```