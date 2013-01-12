# If the out-file is empy it migh imply an error in the running:
#
# FIXME: Does not prodcure error when it should!
if [ -s $FILE_OUT ] ; then
    echo "ok   Installation completed. You may now run $PROGRAM ;) ";
    cat $FILE_OUT;
else
    echo "$FILE_OUT is empty. Writes the error";
    cat $FILE_ERR;
fi ;
exit;