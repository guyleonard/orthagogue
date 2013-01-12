#!/bin/sh -x
echo "Processes file " $1;
cat $1 | awk '
BEGIN {
# initialize all arrays used in for loop
# FIXME: Hvordan oppdatere variabel her??
#    taxon_index=$taxon_index;
    taxon_index=1; 
    count[""]=0;
    fsize[""]=0;
    names[""]=0;
    taxa[""]=0;
}
{
    lines++;
    names[$1]++;
#    print     names[$1] "--" $1;
    if ((x=index($1,"_")) > 0) {
        first = substr($1,1,x-1);
        second = substr($1,x+1,length($1));
# the above is the same as
#        hostname = substr($1,x+1);
        if(first > 0 && second > 0) {
            if(taxon_index == 0) {
                printf("taxon = %s, protein = %s\n", first, second);
                taxa[first]++;
            } else {
                printf("taxon = %s, protein = %s\n", second, first);
                taxa[second]++;
            }
        }
    }
}
END {
    printf("In total %d lines produced.\n", lines);
    printf("The Proteins were:\n");
    for (i in names) {
        printf("%s with %d pairs\n", i, names[i]);
    }
    printf("The Taxa were:\n");
    for (i in taxa) {
        printf("'%s' with %d pairs\n", i, taxa[i]);
    }
}
' 