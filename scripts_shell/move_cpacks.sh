#CURRENT=`pwd`
#BASENAME=`basename $CURRENT`
# Gets the version number by removing the other stuff:
#v_number=`echo $BASENAME | sed s/_[a-zA-Z]*_//`;echo $v_number;
dir_release=$2; #"../$BASENAME.executables";
name_soft=$1;
echo "In 'move_cpacks' dire_release=" $dir_release " name_soft=" $name_soft;
rm -rf $dir_release
mkdir -p $dir_release;
#echo "mkdir -p $dir_release; \n"
mv -f $name_soft-Linux.* _CPack_Packages/ $dir_release; # Moves the executables outside the scope of this package.