#!/usr/bin/perl
# :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Copyright 2012 Ole Kristian Ekseth (oekseth@gmail.com)

# This file is part of heuristical_orthology_building.pl

# heuristical_orthology_building is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# heuristical_orthology_building is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with heuristical_orthology_building. If not, see <http://www.gnu.org/licenses/>.
# :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

#use warnings;
use strict;
use POSIX;
#use types_and_enums qw($format_blastp_evalue $format_blastp_last $format_abc );

# Below a trick forgotten for a c++ programmer (oekseth)
use File::Basename;
use lib basename ($0);
use FindBin;
use lib "$FindBin::Bin";
use File::Copy;
#
use FindBin '$Bin';

# 	if( $? == -1 ){
# 	    print "------------------->>>>>>>>>>>>>>>><<command failed: $!\n";
# 	    exit;
# 	} else {
# 	    printf "command exited with value %d", $? >> 8;
# 	}

#! Return file components:
sub get_file_components {
    my @file_array = split('/', $Bin);    
    if( @file_array > 3 ) {
	my $root = $file_array[@file_array -3 ];
	my @temp = split('_', $root);
	my $software = $temp[0]; my $relase_version = $temp[1];
	my @rel_id = split('\.', $relase_version);
	my $cpack_name = $software . '-';
	my $rel_size = @rel_id;
	for(my $i = 0; $i < $rel_size-1; $i++) {
	    $cpack_name .= $rel_id[$i];
	    if($i < $rel_size -2) {$cpack_name .= '.';}
	}
	$cpack_name .= "-Linux";
# FIXME: remove the below variable!
	print($cpack_name . "\n");
	return ($root, $cpack_name, @rel_id);
    } else {
	return "";
    }
}

sub clean_files_make_cmake {
    my($file_err, @names) = @_;
    foreach my $install_path (@names) {
	chdir($install_path) or die "Can't chdir to $install_path $!"; # In order to have configuration files stored at the same locations.
	my $res2 = `rm -f $file_err *.a *.so CMakeCache.txt Makefile orthAgogue 1>/dev/null 2>/dev/null;`;
	if($install_path ne ".") {chdir("../") or die "Can't chdir to ../ $!";}
    }
}

#
# Builds the directory and (a) moves the files, (b) updates dcoumentation, (c) package the code
sub pre_install {
    my ($root, $cpack_name) = get_file_components();
    #! First call the installation-params to ensure it's the release-version we apply
    system("./install.bash");
    my $dirname = "../" . $root . "_release";
    mkdir($dirname, 0766);
    
    if(1 == 2) {
	system("make package"); # Generates the packages
    } else {
# FIXME: Include below!
	system("make package"); # Generates the packages
	system("cpack -D DEB; mv *.deb ../orthAgogue_release;");
	system("cpack -D RPM; mv *.rpm ../orthAgogue_release;");
	
	my $cpack1 = `cpack -C CPackConfig.cmake;`;
	my $cpack_name_spec = $cpack_name . ".tar.Z";
	move($cpack_name_spec, $dirname);
	$cpack_name_spec = $cpack_name . ".tar.gz";
	move($cpack_name_spec, $dirname);
	$cpack_name_spec = $cpack_name . ".sh";
	move($cpack_name_spec, $dirname);
	$cpack_name_spec = $cpack_name . ".zip";	
	move($cpack_name_spec, $dirname);
	printf("-\tFiles moved to %s\n", $dirname);
}    

#
# Calls some nicely generated scripts:
#printf("at location %s\n", $Bin);
    system("./scripts_shell/remove_temps.sh"); # Remove temporary files.
    system("./scripts_shell/print_status.sh"); # Prints the status.
    system("rm -f report_orthAgogue/*");
    system("mv -f _CPack_Packages/ ../");
    system("mv -f orthAgogue*.tar* ../");
    system("rf -f *.abc *.mci *.map");
    # Remove library files
    system("find . -iname \"*.a\" -exec rm -f {} \;");
    clean_files_make_cmake("err_cleaning.txt", ".");	
    system("tar -cf $dirname/$root.dev.tar ./*");
    system("./scripts_shell/doxy.bash;"); # Produces documentation.
}

sub remove_install_configuration {
    my $FILE_ERR="err_orthAgogue.txt";
    my @names= ("terminal_input", "blast_common", "log_builder", "blast_parsing", "blast_filtering", ".");
    clean_files_make_cmake($FILE_ERR, @names);
}


sub install_at_path {
    my ($install_path, $cmd) = @_;
    chdir($install_path) or die "Can't chdir to $install_path $!"; 
    system($cmd) == 0 or die "Compilation failed: $?";
    my $return_val = system("make");
#    system("make 2>$FILE_ERR");
    if($install_path ne ".") {chdir("../") or die "Can't chdir to ../ $!";}
    return $return_val;
}

#! Iterate the folders had coded in this function, building the script.
sub install_software {
# Sets some properties for the run (easier to use as template then):
    my ($INSTALL_SETTING, $disk_buffer_size, %args_list) = @_;
    my $FILE_OUT="out_orthAgogue.txt";
    my $FILE_ERR="err_orthAgogue.txt";
#    my $cmd2 = "echo \"\" > " . $FILE_OUT; system($cmd2);

    my $PROGRAM="./orthAgogue";
    my ($root, $cpack_name, @version_ids) = get_file_components();
    
    my $v_maJor=$version_ids[0];
    my $v_mInor=$version_ids[1];
    my $v_paTch=$version_ids[2];
    my $log_folder_name="report_orthAgogue";
    my $cmd_add = "";
    foreach my $in (keys %args_list) {
	if($in eq "LOG_FOLDER_NAME") {$log_folder_name = $args_list{"LOG_FOLDER_NAME"};}
	else {
	    $cmd_add = $cmd_add . " " . $in;
	}
    }
    my $ret1 = `module load intelcomp 2>err_module.txt`; # For centOS
    my $ret2 = `module load intelcomp/12.1.0 2>err_module_2.txt`; # For centOS
#    my $rer2=  `make clean`; # For centOS
    my $cmd =  "cmake -D CMAKE_BUILD_TYPE:STRING=$INSTALL_SETTING " 
	. "-D LOG_FOLDER_NAME:STRING=$log_folder_name " 
	. "-D MEMORY_CONSUMPTION_LEVEL=1 "
#	. "-D DISK_BUFFER_SIZE=$disk_buffer_size "
	. "-D VERSION=@version_ids"
	. "-D CPACK_PACKAGE_VERSION_MAJOR=$v_maJor "
	. "-D CPACK_PACKAGE_VERSION_MINOR=$v_mInor "
	. "-D CPACK_PACKAGE_VERSION_PATCH=$v_paTch"
	.  $cmd_add . " .";
    printf("Sends the following line to cmake:\n%s\n\n", $cmd);
    my @names= ( ".");
#    my @names= ("terminal_input", "blast_common", "log_builder", "blast_parsing", "blast_filtering", ".");
    clean_files_make_cmake($FILE_ERR, @names);
    
    # The order of the dependencies the "cmake" file has problems to handle:
#    install_at_path("terminal_input", $cmd);
#    install_at_path("blast_common", $cmd);
    my $res_main = install_at_path(".", $cmd); 
    if(-s $FILE_ERR > 35) {
	printf("--\tDid not complete the run as satsifactory as preferable. (Error-file contains data.) Writes the error for the root (locally found in the file named '%s'):\n", $FILE_ERR);
	printf("--------------------------------------------------------------------------------------\n");
	system("cat $FILE_ERR");
	printf("--------------------------------------------------------------------------------------\n\n");
	printf("--\tIf questions, please contact the developer at oekseth\@gmail.com\n");
    } else {	
	printf("ok   Installation completed. You may now run $PROGRAM ;)\n");
    }
}
my $numArgs = $#ARGV + 1;
if($numArgs == 0) {
    printf("----------------------------------------------------------------------\n");
    printf("-->\tThis message was seen because %d arguments were used. In order to run the software, provide 2 arguments, as described in the line above.\n", $numArgs);
    printf("----------------------------------------------------------------------\n");
} else { # Starts the operations: 

    if($ARGV[0] eq "pre_install") {
	pre_install();
    } elsif($ARGV[0] eq "remove_configs") {
	remove_install_configuration();
    } else {
	my $disk_buffer_size = 0;
	my %args_list = ();
	if($numArgs >= 2) {
	    if(isdigit($ARGV[1])) {$disk_buffer_size = $ARGV[1];} # Sets the disk bufer size if given.
	    if($numArgs >= 3) {
		for(my $i = 2; $i < $numArgs; $i++) {
		# TODO: When complexity increases, improve this part as it comes.
		if($ARGV[$i] eq "-O") {
		    $args_list{"LOG_FOLDER_NAME"} = $ARGV[$i+1];
		} else {
		    $args_list{$ARGV[$i]} = "";
		}
	    }
	    }
	}
	if(($ARGV[0] eq "o2") || ($ARGV[0] eq "RELEASE") || ($ARGV[0] eq "o3") || ($ARGV[0] eq "fart") || ($ARGV[0] eq "fast")) {
	    my $INSTALL_SETTING = "RELEASE";
	    if(!$disk_buffer_size) {$disk_buffer_size=1024*1024*20;} # Sets the disk buffer size:
	    
	    install_software($INSTALL_SETTING, $disk_buffer_size, %args_list);
	} elsif(($ARGV[0] eq "DEBUG") || ($ARGV[0] eq "bug")) {
	    my $INSTALL_SETTING = "DEBUG";
	    if(!$disk_buffer_size){$disk_buffer_size=2000;} # Sets the disk buffer size:
	    install_software($INSTALL_SETTING, $disk_buffer_size, %args_list);
	} else {
	    printf("!!\tAn error, as the input given='%s' did not correspond to the parameters set.\n", $ARGV[0]);
	}
    }
} 

1;
