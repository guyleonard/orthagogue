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

use warnings;
use strict;
use POSIX;
#use types_and_enums qw($format_blastp_evalue $format_blastp_last $format_abc );

# Below a trick forgotten for a c++ programmer (oekseth)
use File::Basename;
use lib basename ($0);
use FindBin;
use lib "$FindBin::Bin";
use File::Copy;

# The local modules:
use blast_parsing;
use blast_filtering;
use get_orthologs;
use get_inparalogs;
use get_co_orthologs;
# Defines the booleans in the code, making it easier to shift this code to c++ (as always is or path of walk).
use constant false => 0;
use constant true  => 1;
# Defines the booleans in the code, making it easier to shift this code to c++ (as always is or path of walk).
use constant false => 0;
use constant true  => 1;
# # Defines the file format to parse:
use constant format_blastp_evalue   => 'blastp_c_11';
use constant format_blastp_last     => 'blastp_c_12';
use constant format_abc             => 'abc';


my $control_output_folder = "control_set/";

my $USE_MPI = false;

sub control_files {
    my($blastp_file, $omcl_path, $turbo_path) = @_;
#
# The parsing of the blastp file used for both omcl and to:
    my %taxa = blast_parsing::get_blastp_info($blastp_file, format_blastp_evalue);
    my %blastp= blast_parsing::get_hash($blastp_file, format_blastp_evalue, #only_best_alignment=
					1); # Note: this value is of high imporatance!
#				    0); # Note: this value is of high imporatance!
#    blast_parsing::print_taxa(\%taxa);

    # Gets the raw data in order to visually verify the process:
    my %blastp_raw = blast_parsing::get_raw_blastp_file($blastp_file, format_blastp_evalue);
    my %blastp_average = blast_parsing::make_avarage(\%blastp, \%blastp_raw);
    if(!(defined \%blastp_raw)) { # Included to avoid garbage collection, wich sometimes happends.
	printf("An error in the code: blast_raw hash not defined!\n");
    }
    blast_filtering::write_averaged_data_file($control_output_folder, \%blastp_average);

#
# The ortholog procedure:
    my %inparalog_limit = ();
    my %orthologs = get_orthologs::build_orthologs(\%blastp_average, \%inparalog_limit);
    my $omcl_file_orthologs = $omcl_path . "orthologs.abc.sort";
    my $turbo_file_orthologs = $turbo_path . "orthologs.abc";
    my %omcl = ();
    if($omcl_file_orthologs ne "") {%omcl = blast_parsing::get_hash($omcl_file_orthologs, format_abc);}
    my %turbo= ();
    if($turbo_file_orthologs ne "") {%turbo = blast_parsing::get_hash($turbo_file_orthologs, format_abc);}
    get_orthologs::write_ortholog_control_file($control_output_folder, \%orthologs);
    get_orthologs::compare_orthologs($control_output_folder, \%orthologs, \%omcl, \%turbo, \%blastp_average, \%blastp_raw);

#    
# The inparalog procedure:
    if(!(defined \%inparalog_limit)) { # Included to avoid garbage collection, wich sometimes happends.
	printf("An error in the code: inparalog_limit hash not defined!\n");
    }
    my %inparalogs = get_inparalogs::build_inparalogs(\%blastp_average, \%inparalog_limit);
    my $omcl_file_inparalogs = $omcl_path . "inparalogs.abc.sort";
    my $turbo_file_inparalogs = $turbo_path . "inparalogs.abc";
    my %omcl_inpa = ();
    if($omcl_file_inparalogs ne "") {%omcl_inpa = blast_parsing::get_hash($omcl_file_inparalogs, format_abc);}
    my %turbo_inpa= ();
    if($turbo_file_inparalogs ne "") {%turbo_inpa = blast_parsing::get_hash($turbo_file_inparalogs, format_abc);}
    get_inparalogs::write_inparalog_control_file($control_output_folder, \%inparalogs);
    get_inparalogs::compare_inparalogs($control_output_folder, \%inparalogs, \%omcl_inpa, \%turbo_inpa, \%blastp_average, \%blastp_raw, \%inparalog_limit);

#    
# The co-ortholog procedure:
    my %co_orthologs = get_co_orthologs::build_co_orthologs(\%blastp_average, \%inparalogs, \%orthologs,
							    1 # The nss: must be set to '1' in order to be able tom compare with OrthoMcl
							    );
    my $omcl_file_co_orthologs = $omcl_path . "coorthologs.abc.sort";
    my $turbo_file_co_orthologs = $turbo_path . "co_orthologs.abc";
    my %omcl_co_orthologs = ();
    if($omcl_file_co_orthologs ne "") {%omcl_co_orthologs = blast_parsing::get_hash($omcl_file_co_orthologs, format_abc);}
    my %turbo_co_orthologs= ();
    if($turbo_file_co_orthologs ne "") {%turbo_co_orthologs = blast_parsing::get_hash($turbo_file_co_orthologs, format_abc);}
    get_co_orthologs::write_co_ortholog_control_file($control_output_folder, \%co_orthologs);
    get_co_orthologs::compare_co_orthologs($control_output_folder, \%co_orthologs, \%omcl_co_orthologs, \%turbo_co_orthologs, \%blastp_average, \%blastp_raw);
}

#! Test the settings: 
sub test_settings {
    my $validate_folder = "validated_output";
    my $folder_exsists = blast_filtering::make_folder_if_not_set($validate_folder);
    my %tests; 
    if($USE_MPI == false) {
	%tests = (
		  simple =>  "./orthAgogue -i empirical_tests/goodProteins.blast"
 ,
 		  simple_dbs_1024 =>  "./orthAgogue -i empirical_tests/goodProteins.blast -dbs 1024",
 		  simple_dbs_10024 =>  "./orthAgogue -i empirical_tests/goodProteins.blast -dbs 10024",
 		  simple_dbs_40024 =>  "./orthAgogue -i empirical_tests/goodProteins.blast -dbs 40024",
 		  norm_all =>  "./orthAgogue -i empirical_tests/goodProteins.blast -A",
 		  simple_e0 =>  "./orthAgogue -i empirical_tests/goodProteins.blast -e 0",
 		  simple_e40 =>  "./orthAgogue -i empirical_tests/goodProteins.blast -e 40",
		  simple_050 =>  "./orthAgogue -i empirical_tests/goodProteins.blast -o 50",
		  simple_090 =>  "./orthAgogue -i empirical_tests/goodProteins.blast -o 90",
 		  simple_nii =>  "./orthAgogue -i empirical_tests/goodProteins.blast -sco",
 		  simple_nss =>  "./orthAgogue -i empirical_tests/goodProteins.blast -bho",
 		  # simple_nn =>  "./orthAgogue -i empirical_tests/goodProteins.blast -nn",
 		  # simple_nn_c2 =>  "./orthAgogue -i empirical_tests/goodProteins.blast -nn -c 2",
 		  # simple_nn_c3 =>  "./orthAgogue -i empirical_tests/goodProteins.blast -nn -c 3",
 		  # simple_nn_c4 =>  "./orthAgogue -i empirical_tests/goodProteins.blast -nn -c 4",
 		  # simple_nn_c5 =>  "./orthAgogue -i empirical_tests/goodProteins.blast -nn -c 5"
	    
		  );
    } else {
	%tests = (
 		  simple =>  "mpirun -np 2 ./orthAgogue -i empirical_tests/goodProteins.blast -dbs 8000",
  		  simple_dbs_1024 =>  "mpirun -np 2 ./orthAgogue -i empirical_tests/goodProteins.blast -dbs 1024",
  		  simple_dbs_10024 =>  "mpirun -np 2 ./orthAgogue -i empirical_tests/goodProteins.blast -dbs 10024",
  		  simple_dbs_40024 =>  "mpirun -np 2 ./orthAgogue -i empirical_tests/goodProteins.blast -dbs 10024",
  		  norm_all =>  "mpirun -np 2 ./orthAgogue -i empirical_tests/goodProteins.blast -A -dbs 10024",
  		  simple_e0 =>  "mpirun -np 2 ./orthAgogue -i empirical_tests/goodProteins.blast -e 0 -dbs 10024",
  		  simple_e40 =>   "mpirun -np 2 ./orthAgogue -i empirical_tests/goodProteins.blast -e 40 -dbs 10024",
  		  simple_050 =>   "mpirun -np 2 ./orthAgogue -i empirical_tests/goodProteins.blast -o 50 -dbs 10024",
  		  simple_090 =>   "mpirun -np 2 ./orthAgogue -i empirical_tests/goodProteins.blast -o 90 -dbs 10024 ",
  		  simple_nii =>   "mpirun -np 2 ./orthAgogue -i empirical_tests/goodProteins.blast -sco -dbs 10024 ",
  		  simple_nss =>   "mpirun -np 2 ./orthAgogue -i empirical_tests/goodProteins.blast -bho -dbs 10024 ",
  		  simple_nn =>    "mpirun -np 2 ./orthAgogue -i empirical_tests/goodProteins.blast -nn -dbs 10024 ",
   		  simple_nn_c2 => "mpirun -np 2 ./orthAgogue -i empirical_tests/goodProteins.blast -nn -c 2 -dbs 10024 ",
   		  simple_nn_c3 => "mpirun -np 2 ./orthAgogue -i empirical_tests/goodProteins.blast -nn -c 3 -dbs 10024 ",
   		  simple_nn_c4 => "mpirun -np 2 ./orthAgogue -i empirical_tests/goodProteins.blast -nn -c 4 -dbs 10024 ",
  		  simple_nn_c5 => "mpirun -np 4 ./orthAgogue -i empirical_tests/goodProteins.blast -nn -c 5 -dbs 1024 ",		 

  		  small_11M_nn_c1_b =>  "mpirun -np 1 ./orthAgogue -i empirical_tests/file_11m.blast -nn -c 1 -dbs 50024 -t 1 -p 0 -s '_'", 
  		  small_11M_nn_c1_n1 => "mpirun -np 1 ./orthAgogue -i empirical_tests/file_11m.blast -nn -c 1 -dbs 50024 -t 1 -p 0 -s '_'", 
  		  small_11M_nn_c1_n5 => "mpirun -np 5 ./orthAgogue -i empirical_tests/file_11m.blast -nn -c 1 -dbs 50024 -t 1 -p 0 -s '_'", 
 		  small_11M_nn_c1 => "mpirun -np 3 ./orthAgogue -i empirical_tests/file_11m.blast -nn -c 5 -dbs 50024 -t 1 -p 0 -s '_'" 

		  );
    }
    my $test_1_cmd = "./orthAgogue -i empirical_tests/goodProteins.blast";
    my $test_1_name = "test_1";
    foreach my $in (keys %tests) {
	my $test_path = $validate_folder . "/" . $in;
	system($tests{$in}); # Builds the file to compare
	my $retval = open(my $file_res, "<", $test_path);
	if(defined $retval) { # Template exsists: compare:
	    close($file_res) or die $!;
	    my %correct = blast_parsing::get_hash($test_path, format_abc, #only_best_alignment=
						  0); # Note: this value is of high imporatance!	    
	    my %new_f = blast_parsing::get_hash("all.abc", format_abc, #only_best_alignment=
						0); # Note: this value is of high imporatance!	    
	    my $cnt_differences = blast_filtering::compare_verbouse($in, \%correct, \%new_f);
	    if($cnt_differences == 0) {
#		printf("success:\t%s\n", $tests{$in});
		use Term::ANSIColor;
#		print colored ("Text\n", 'bold blue on_white');  
		print colored ("success:\t", 'dark bold green');
#		print "success:\t";
#		print color 'reset';
		printf("%s\n", $tests{$in});

	    } else {
		printf("!!\tAborts due to differences: Please contact the developer at oekseth\@gmail.com if this error is seen: %s \n",  $tests{$in});
		printf("PS: If this error was during the development pahse, remember also validating the template named %s\n", $test_path);
		exit(2);
	    }
	} else { # Does not exsists: Makes the generated file our new template
	    printf("\tMakes a new template named %s\n", $test_path);
	    move("all.abc", $test_path) or die $!;
	}
    }
}
# FIXME: validate the below line!!
#./../orthAgogue -i goodProteins.blast -nn -O turbo/gp_file/ -pd -bho 1>results_detailed.txt 2>err.txt; 
#my $correct_turbo_params = `./../orthAgogue -i goodProteins.blast -nn -O turbo/ -pd -bho 1>results_detailed.txt 2>err.txt`;
#$correct_turbo_params = `./../orthAgogue -i goodProteins.blast -nn -O turbo/gp_file/ -pd`;

#my %blastp= get_hash($blastp_file, format_blastp_last);
#my $blastp_file = "sample_test_for_max_filtering.blast";
my $omcl_file = "";
my $blastp_file = "goodProteins.blast";
my $turbo_file = "turbo/gp_file/"; #orthologs.abc";
$omcl_file = "omcl/gp_file/"; #orthologs.abc.sort";

#my $turbo_file = "turbo/gp_file/orthologs.abc";
#my $blastp_file = "all_vm.blast";
#my $turbo_file = "turbo/all_vm/";
#my $omcl_file = "omcl/gp_file/orthologs.abc.sort";
#my $turbo_file = "turbo/gp_file/orthologs.abc";
#control_files($blastp_file, $omcl_file, $turbo_file);
my $numArgs = $#ARGV + 1;
if($numArgs !=3 && $numArgs != 2 && $numArgs != 4) {
    printf("----------------------------------------------------------------------\n");
    printf("Building a control set containing filtered relations, storing detailed data in folder: %s\n", $control_output_folder);
    printf("-\tValidates implementations found in softwares, eg, OrthoMcl and OrthaGogue\n");
    printf("-\tUsage: %s <blastp_file> <OrthaGogue> <optionally argument; other software implementations, ie, OrthoMcl>\n", $0);
    printf("-->\tThis message was seen because %d arguments were used. In order to run the software, provide either 2 or 3 arguments, as described in the line above.\n", $numArgs);
    printf("\nThe software was developed by O.K. Ekseth under supervison of Dr. V.N. Mironov. Questions to be forwarded to oekseth\@gmail.com\n");
    printf("----------------------------------------------------------------------\n");
} else { # Starts the operations: 
#    foreach my $argnum (0 .. $#ARGV) {print "[$argnum]\t$ARGV[$argnum]\n";}
    if($numArgs == 2) {control_files($ARGV[0], "", $ARGV[1]);} #$omcl_file, $turbo_file);   
    elsif($numArgs == 3 || $numArgs == 4) {
	if($numArgs == 4) {
	    if($ARGV[3] eq "MPI") {
		$USE_MPI = true;
		printf("Uses MPI\n");
	    } 
	}
	control_files($ARGV[0], $ARGV[2], $ARGV[1]);
    } #$omcl_file, $turbo_file);   
    test_settings();
} 

#control_files($blastp_file, $omcl_file, $turbo_file);

1; # The return value;
