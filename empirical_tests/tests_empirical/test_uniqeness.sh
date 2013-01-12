#!/bin/sh
echo ""
make; ./orthaGogue -i all_vm.blast -c 1; # Generates the output
echo "Tests uniqueness for the co-orthologs:"
cut -f 1,2 co_orthologs.abc | sort > co_sorted.abc; # Sorting allowing duplicates.
cut -f 1,2 co_orthologs.abc | sort -u > co_unique_sorted.abc; # Sorting, removing duplicates.
diff co_sorted.abc co_unique_sorted.abc; # If successfull, no output is genereated.
wc -l co_sorted.abc co_unique_sorted.abc; # If successfull, resulting numbers are equal.
echo ""
echo "Tests uniqueness for the ortholog file:"
cut -f 1,2 orthologs.abc | sort > orthologs_sorted.abc; # Sorting orthologs, alowing duplicates.
cut -f 1,2 orthologs.abc | sort -u > orthologs_unique_sorted.abc; # Sorting, removing duplicates.
diff orthologs_sorted.abc orthologs_unique_sorted.abc; # If successfull, no output is generated.
wc -l orthologs_sorted.abc orthologs_unique_sorted.abc; # If successfull, resulting numbers are equal.
echo ""
echo "Tests uniqueness for the complete output file:"
cut -f 1,2 all.abc | sort > all_sorted.abc; # Sorting allowing duplicates.
cut -f 1,2 all.abc | sort -u > all_unique_sorted.abc; # Sorting, removing duplicates.
diff all_sorted.abc all_unique_sorted.abc; # If successfull, no output is generated.
wc -l all_sorted.abc all_unique_sorted.abc; # If successfull, resulting numbers are equal.
exit;