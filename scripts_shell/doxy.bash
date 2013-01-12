#!/bin/bash

#------------------------------------------------------------
# GENERATES THE DOCUMENTATION
#------------------------------------------------------------
CURRENT=`pwd`
BASENAME=`basename $CURRENT`
# Gets the version number by removing the other stuff:
v_number=`echo $BASENAME | sed s/_[a-zA-Z]*_//`;
#doc_path="..\/..\/Documentation";  # Backtics to avoid misunderstand in the 'sed' interpretator
doc_path="..\/Documentation";  # Backtics to avoid misunderstand in the 'sed' interpretator
#doc_path="../../Documentation";  
#echo "Buidls doxy with v_number=" $v_number;

#echo "PWD=" $CURRENT;
# Replaces the varaibel included in the Doxyfile with the version number found above:
sed s/project_number/$v_number/ < Doxyfile > temp1_Doxyfile;
sed s/doc_path/$doc_path/ < temp1_Doxyfile > temp_Doxyfile;
#sed s|doc_path|$doc_path| < temp1_Doxyfile > temp_Doxyfile;
#echo "sed s/doc_path/$doc_path/ < temp1_Doxyfile > temp_Doxyfile;"
echo "starts the building of the doxy.."
doxygen temp_Doxyfile 1>out_doxy.txt 2>err_doxy.txt;cat err_doxy.txt;
resu=`rm temp_Doxyfile temp1_Doxyfile`;
# Remove the temporary files:
#rm temp_Doxyfile out_doxy.txt temp1_Doxyfile;
PATH="../Documentation/latex/"; 
if [ -d "$PATH" ]; then
    cd $PATH;
    # Control will enter here if $DIRECTORY exists
else
    cd "Documentation/latex/";
fi
CURRENT_TEMP=`pwd`;
echo "pdflatex in $CURRENT_TEMP";
#result=`pdflatex refman.tex 1>out_latex.txt`;  #2>err_latex.txt 
result=`/usr/bin/pdflatex refman.tex 1>out_latex.txt`;  #2>err_latex.txt 
result=`/usr/bin/pdflatex refman.tex 1>out_latex.txt`;  #2>err_latex.txt 
#cat err_latex.txt;
echo $result;
cd $CURRENT;
exit;