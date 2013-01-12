#!/bin/sh
##!/bin/sh  -x
if [ "$#" -eq 1 ] ; then 
    INPUT_FILE=$1;
#    INPUT_FILE="all_vm.blast"; 
    cmd="./orthaGogue -i $INPUT_FILE --no_normalization --not_sum_scores --use_last_column --output_dir .";
    echo $cmd;
#`$cmd`;
    ABC_FILE="all.abc"
    echo "'''''''''''''''''''''''''''''''''''''''''''''''''''''''"
    echo "Script written for file $INPUT_FILE $ABC_FILE to verify the \"--not_sum_scores\" option in TurboOrtho."
    echo "-------------------------------------------------------"
    cmd="./tests_empirical/print_reciprocal_pair.sh \$one \$two $INPUT_FILE $ABC_FILE;";
    cmd_capsulated="./tests_empirical/print_reciprocal_pair.sh";
    # Note: In order to use global variables, <'> is used to quote the global (dynamic) variables:
    echo "#!/bin/sh" > temp.sh;
    awk '{print "'$cmd_capsulated' \"" $1 "\" \"" $2 "\" '$INPUT_FILE' '$ABC_FILE'"}' $ABC_FILE >>temp.sh
    chmod 711 temp.sh; ./temp.sh; rm temp.sh;
else
    echo "!!  Not able to validation as blastp file not specified (requires 1 input, but  " $# " arguments included)"
fi

exit;