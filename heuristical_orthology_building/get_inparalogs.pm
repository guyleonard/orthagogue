#!/usr/local/bin/perl
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

package get_inparalogs;
use warnings;
use strict;
use POSIX;
#use types_and_enums qw($format_blastp_evalue $format_blastp_last $format_abc );
use blast_parsing;
use blast_filtering;
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

# Names of the files generated:
my $inparalog_name = "inparalogs_control.abc";

# The output names of the differences:
my $name_fn_gogue = "inparalogs_FN_gogue";
my $name_fp_gogue = "inparalogs_FP_gogue";
my $name_fn_omcl = "inparalogs_FN_omcl";
my $name_fp_omcl = "inparalogs_FP_omcl";

# Settings changing default output:
my $print_disregarded_inparalogs = 0; # If set to '1' prints the possibile inparalogs discraded.
my $print_inparalog_limits = 0;

#! Builds the orthologs
sub build_inparalogs {
    my ($blastp_, $inparalog_limit_) = @_;
    my %blastp = %$blastp_; my %inparalog_limit = %$inparalog_limit_;
    my %inparalogs;
    # Starts building a non-reciprocal ortholog set:    
    foreach my $in (keys %blastp) {
	if(defined($inparalog_limit{$in})) {
	    my $taxon_in = blast_parsing::get_taxon_type($in);
	    foreach my $taxon_out (keys %{ $blastp{$in} }) {
		# Only build relations for those combinations not constituting inparalogs:
		if($taxon_in eq $taxon_out) {
		    while (my ($out, $value) = each %{ $blastp{$in}{$taxon_out} } ) {		    		    
			if($in ne $out) { # They are not a self-comparison: 
			    if(defined($inparalog_limit{$out})) {
				# Below procedure ensures reciprocability:
				if(($value >= $inparalog_limit{$in}) && ($value >= $inparalog_limit{$out})) {
				    # NOTE: The use of 'taxon_id' as variable, enable using of default proecedure, even though it really has no point, if you ignore the scope of readability.
				    $inparalogs{$in}{$taxon_out}{$out} = $value;
				    $inparalogs{$in}{$taxon_out}{$out} = $value;
				    if($print_disregarded_inparalogs) {printf("(disregarded)\t(%s %s), with value %f < max(%f, %f)\n",$in, $out, $value,  $inparalog_limit{$in},  $inparalog_limit{$out});}
				}
			    }
			}
		    }
		}
	    }
	} 
    }
    return %inparalogs;
}


sub write_inparalog_control_file {
   my ($control_output_folder, $list_) = @_;
   blast_filtering::write_control_file($control_output_folder, $inparalog_name, $list_);
}

sub print_preferences {
    my ($id_2_param, $outsNotFoundInList_, $list_, $blastp_average_, $blastp_raw_, $inparalog_limit_) = @_;
    if(!(defined $blastp_raw_)) {
	printf("Not defined the blastp_raw; an error\n");
    }
    if(!(defined $blastp_average_)) {
	printf("Not defined the blastp_average_; an error\n");
    }
#    my %blastp_average = %$blastp_average_;
#    %blastp_raw = %$blastp_raw_;
    # Remaps:
    my %outsNotFoundInList = %$outsNotFoundInList_; my %list = %$list_; my %blastp_average = %$blastp_average_;
    my %inparalog_limit = %$inparalog_limit_;
    foreach my $in (keys %outsNotFoundInList) {
	my $taxon_in = blast_parsing::get_taxon_type($in);
	foreach my $taxon_out (keys %{ $outsNotFoundInList{$in} }) {
	    my $size = scalar(keys %{ $list{$in}{$taxon_out} });
	    if($size > 0) {
		printf("(**)\t Inparalogs for %s:\n", $in);
		while (my ($out, $value) = each %{ $outsNotFoundInList{$in}{$taxon_out} } ) {
		    blast_parsing::print_raw_data_recip($blastp_raw_, $in, $out, $taxon_in, $taxon_out, $blastp_average_, $inparalog_limit{$in},  $inparalog_limit{$out});
		  }
		printf("--> \t But instead %s prefered:\n", $id_2_param);
		while (my ($out, $value) = each %{ $list{$in}{$taxon_out} } ) {
		    blast_parsing::print_raw_data_recip($blastp_raw_, $in, $out, $taxon_in, $taxon_out, $blastp_average_, $inparalog_limit{$in},  $inparalog_limit{$out});
		  }
		printf("-----------------------------------------\n");
	    } else {
		printf("\t(NULL) The inparalogs were not defined in %s, but should have been:\n", $id_2_param); 
		while (my ($out, $value) = each %{ $outsNotFoundInList{$in}{$taxon_out} } ) {
		    blast_parsing::print_raw_data_recip($blastp_raw_, $in, $out, $taxon_in, $taxon_out, $blastp_average_, $inparalog_limit{$in},  $inparalog_limit{$out});
		  }
		printf("-----------------------------------------\n");
	    }
	} 
    }
}

#! Compares OrthoMcl and TurboOrtho to the control set:
sub compare_control_inparalogs_to_other {
    my($other_id, $control_output_folder, $inparalogs_, $other_inparalogs_, $blastp_average_, $name_fp_other, $name_fn_other, $blastp_raw_, $inparalog_limit_) = @_;
    if(!(defined $blastp_raw_)) {
	printf("Not defined the blastp_raw (1); an error\n");
    }
    # The returned list holds those inparalogs not found in orthaGogue's dataset:
    my($TP_other_1, $FN_other, %inparalogs_fn_other) = blast_filtering::get_difference($inparalogs_, $other_inparalogs_);
    # Writes the lacking inparalogs to an own file:
    blast_filtering::write_control_file($control_output_folder, $name_fn_other, \%inparalogs_fn_other);

    # The returned list holds those inparalogs not found in controls dataset:
    my($TP_other_2, $FP_other, %inparalogs_fp_other) = blast_filtering::get_difference($other_inparalogs_, $inparalogs_);
    # Writes those errate inparalogs to an own file:
    blast_filtering::write_control_file($control_output_folder, $name_fp_other, \%inparalogs_fp_other);

    if($FN_other > 0) {
	printf("- The preferences are for the %d FN in %s are:\n", $FN_other, $other_id);	
	print_preferences($other_id, \%inparalogs_fn_other, $other_inparalogs_, $blastp_average_, $blastp_raw_, $inparalog_limit_);
    }
#printf("- The preferences are for the %d false hits (FP) in %s are:\n", $FP_other, $other_id);	
    if($FP_other > 0) {
	printf("- The preferences are for the %d false hits (FP) in %s are:\n", $FP_other, $other_id);	
	print_preferences("control", \%inparalogs_fp_other, $inparalogs_, $blastp_average_, $blastp_raw_, $inparalog_limit_);
    }
    if($TP_other_1 != $TP_other_2) {
	printf("- For %s, we observe that wile searsing for the FP we find %d TP or %d TP if searching for FP.\n", $other_id, $TP_other_1, $TP_other_2);
    }
    return ($TP_other_1, $FP_other, $FN_other);
}
#! Compares OrthoMcl and TurboOrtho to the control set:
sub compare_inparalogs {
    my($control_output_folder, $inparalogs_, $omcl_inparalogs_, $gogue_inparalogs_, $blastp_average_, $blastp_raw_, $inparalog_limit_) = @_;
    if(!(defined $blastp_raw_)) {	printf("the blastp_raw is not defined: an error in the file.\n");}
    if($print_inparalog_limits) {blast_filtering::print_inparalog_limits($inparalog_limit_);}
# Tests TurboOrtho:
    my($TP_gogue, $fp_gogue, $fn_gogue) = compare_control_inparalogs_to_other("gogue", $control_output_folder, $inparalogs_, $gogue_inparalogs_, $blastp_average_, $name_fp_gogue, $name_fn_gogue, $blastp_raw_, $inparalog_limit_);
# Tests OrthoMcl:
    my($TP_omcl, $fp_omcl, $fn_omcl) = (0,0, 0);
    if((%$omcl_inparalogs_)) {
	($TP_omcl, $fp_omcl, $fn_omcl)    = compare_control_inparalogs_to_other("omcl", $control_output_folder, $inparalogs_, $omcl_inparalogs_, $blastp_average_, $name_fp_omcl, $name_fn_omcl, $blastp_raw_, $inparalog_limit_);
    } else {printf("...\t(Does not make a control of the omcl results for the inparalogs.)\n");}
    printf("\n#\tInparalogs: In sum we see that \n-\torthaGogue has %d TP, %d FP and %d FN, while\n-\tOrthoMcl has %d TP, %d FP and %d FN.\n", $TP_gogue, $fp_gogue, $fn_gogue, $TP_omcl, $fp_omcl, $fn_omcl);
}


1;
