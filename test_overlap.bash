#!/bin/bash
#
# rm *.abc *.mci; 
# make 1>out.txt; ./orthAgogue -i empirical_tests/goodProteins.blast --print_blast_filter_data; # > out_without_overlap_thrs.txt
# wc -l *.abc;
# rm *.abc *.mci; 
# make 1>out.txt; ./orthAgogue -i empirical_tests/goodProteins.blast --print_blast_filter_data -lc; # > out_without_overlap_thrs.txt
# wc -l *.abc;

rm *.abc *.mci; 
echo "./orthAgogue -i /norstore/project/ssb/workspace/data/blast/cco.blast -p 0 -t 1 -s '_' --print_blast_filter_data;\n";
./orthAgogue -i /norstore/project/ssb/workspace/data/blast/cco.blast -p 0 -t 1 -s '_' --print_blast_filter_data; # > out_without_overlap_thrs.txt
wc -l *.abc; rm *.abc *.mci; 
echo "./orthAgogue -i /norstore/project/ssb/workspace/data/blast/cco.blast -p 0 -t 1 -s '_' -o 50 --print_blast_filter_data;\n";
./orthAgogue -i /norstore/project/ssb/workspace/data/blast/cco.blast -p 0 -t 1 -s '_' -o 50 --print_blast_filter_data; # > out_o50_overlap_thrs.txt
wc -l *.abc; rm *.abc *.mci; 
echo "./orthAgogue -i /norstore/project/ssb/workspace/data/blast/cco.blast -p 0 -t 1 -s '_' -o 90 --print_blast_filter_data;\n";
./orthAgogue -i /norstore/project/ssb/workspace/data/blast/cco.blast -p 0 -t 1 -s '_' -o 90 --print_blast_filter_data; # > out_o90_overlap_thrs.txt
wc -l *.abc; rm *.abc *.mci; 
echo "./orthAgogue -i /norstore/project/ssb/workspace/data/blast/cco.blast -p 0 -t 1 -s '_' -lc --print_blast_filter_data;\n";
./orthAgogue -i /norstore/project/ssb/workspace/data/blast/cco.blast -p 0 -t 1 -s '_' -lc --print_blast_filter_data; 
wc -l *.abc; rm *.abc *.mci; 
echo "./orthAgogue -i /norstore/project/ssb/workspace/data/blast/cco.blast -p 0 -t 1 -s '_' -o 50 -lc --print_blast_filter_data;\n";
./orthAgogue -i /norstore/project/ssb/workspace/data/blast/cco.blast -p 0 -t 1 -s '_' -o 50 -lc --print_blast_filter_data; 
wc -l *.abc; rm *.abc *.mci; 
echo "./orthAgogue -i /norstore/project/ssb/workspace/data/blast/cco.blast -p 0 -t 1 -s '_' -o 90 -lc --print_blast_filter_data;\n";
./orthAgogue -i /norstore/project/ssb/workspace/data/blast/cco.blast -p 0 -t 1 -s '_' -o 90 -lc --print_blast_filter_data; 
wc -l *.abc;
