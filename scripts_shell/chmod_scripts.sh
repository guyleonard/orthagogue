#!/bin/sh
ROOT=".";
# Use the following commands to change file or folder permissions:
# chmod (change file modes)
# chown (change file owner)
# chgrp (change file group owner)

# The following letters represent
# " u " - user/owner
# " g " - group owner
# " o " - all other users
# " a " - for all: user/owner, group owner and all other users
# " r " - read permission
# " w " - write permission
# " x " - execute permission 

# Chmods the scripts in the given folder, thereby avoid using the "sh <file_name>.sh" command.
for scripts in $(find $ROOT -name "*.sh") ; do
    chmod -f ugo+x $scripts;
done

exit;