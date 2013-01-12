#!/bin/sh 
file1="omcl/coorthologs.abc.sort";
awk '
BEGIN {
file2="turbo/orthologs.abc";
#  file="goodProteins.blast";
#  system("head " file)
#  system("head goodProteins.blast");
#  found=system("grep \"At1g01190_at\" $file2");
  echo $found;
}
{
#  printf("%s %s %f\n", $1, $2, $3); 
#  found=system("grep \"At1g01190_at\" " file2);
 # printf("found=%s\n", $found);
  found=system("grep \"" 1 "\" "file2);
#  echo $found;
#  result=`grep $0`;
} END {
printf("At the end\n");
} 
' $file1;
#grep "$1" $3 | grep "$2" | sort --key=1,2 -u | awk -F " " '{ print "(blastp)  " $1 "\t" $2 "\t" $12}';# #> $fname;
#grep "$1" $4 | grep "$2"|  awk -F " " '{ print "(turbOO)  " $1 "\t" $2 "\t" $3}';# #> $fname;
#echo "-------------------------------------------------------"
