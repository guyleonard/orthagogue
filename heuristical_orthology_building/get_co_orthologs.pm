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

package get_co_orthologs;
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
my $co_ortholog_name = "co_orthologs_control.abc";

# The output names of the differences:
my $name_fn_gogue = "co_orthologs_FN_gogue";
my $name_fp_gogue = "co_orthologs_FP_gogue";
my $name_fn_omcl = "co_orthologs_FN_omcl";
my $name_fp_omcl = "co_orthologs_FP_omcl";

# Settings changing default output:
my $print_disregarded_co_orthologs = 0; # If set to '1' prints the possibile co_orthologs discraded.
#my $print_inparalog_limits = 0;

# Adds inparalogs for a given set of orthologs.
sub add_co_orthologs_inparalogs_inparalogs {
    my ($blastp_, $inparalogs_, $orthologs_, $co_orthologs, $in, $taxon_in, $out, $taxon_out) = @_;
    my %blastp = %$blastp_; my %inparalogs = %$inparalogs_; my %orthologs = %$orthologs_;
    # Adds the the inparalogs of 'in' to the 'out' relation of the ortholog pair:
#    my $size = scalar(keys %{$inparalogs{$in}{$taxon_in } });
#    printf("starts 'add_co_orthologs_simple(..) with size %d and in=%s, taxon_in=%s\n", $size, $in, $taxon_in);
    foreach my $inpa_left (keys %{ $inparalogs{$in}{$taxon_in}}) {
	foreach my $inpa_right (keys %{ $inparalogs{$out}{$taxon_out}}) {
	    
	    if(defined($blastp{$inpa_right}{$taxon_in}{$inpa_left})) {
		if(!(defined $orthologs{$inpa_right}{$taxon_in}{$inpa_left})) {
		    # A brain-dead way of ensuring uniqeness:
		    if(!(defined $co_orthologs->{$inpa_right}{$taxon_in}{$inpa_left})) {
			$co_orthologs->{$inpa_right}{$taxon_in}{$inpa_left} =  $blastp{$inpa_right}{$taxon_in}{$inpa_left};
			$co_orthologs->{$inpa_left}{$taxon_out}{$inpa_right} = $blastp{$inpa_right}{$taxon_in}{$inpa_left};
		    }
#		printf("(inserted)\t %s %s \n", $out, $inpa); printf("(inserted)\t %s %s \n", $inpa, $out);
		} else {
		    if($print_disregarded_co_orthologs) {
			printf("(discarded_due_to_ortholog)\t %s %s \n", $inpa_left, $inpa_right);
			printf("(discarded_due_to_ortholog)\t %s %s \n", $inpa_right, $inpa_left);
		    }
		}
	    } else {
		if($print_disregarded_co_orthologs) {
		    printf("(discarded_due_to_not_found)\t %s %s \n", $inpa_left, $inpa_right);
		    printf("(discarded_due_to_not_found)\t %s %s \n", $inpa_right, $inpa_left);
		}
	    }
	}
    }
}

# Adds inparalogs for a given ortholog.
sub add_co_orthologs_inparalogs_ortholog {
    my ($blastp_, $inparalogs_, $orthologs_, $co_orthologs, $in, $taxon_in, $out, $taxon_out) = @_;
    my %blastp = %$blastp_; my %inparalogs = %$inparalogs_; my %orthologs = %$orthologs_;
    # Adds the the inparalogs of 'in' to the 'out' relation of the ortholog pair:
#    my $size = scalar(keys %{$inparalogs{$in}{$taxon_in } });
#    printf("starts 'add_co_orthologs_simple(..) with size %d and in=%s, taxon_in=%s\n", $size, $in, $taxon_in);
    foreach my $inpa (keys %{ $inparalogs{$in}{$taxon_in}}) {
	if(defined($blastp{$out}{$taxon_in}{$inpa})) {
	    if(!(defined $orthologs{$out}{$taxon_in}{$inpa})) {
		$co_orthologs->{$out}{$taxon_in}{$inpa} = $blastp{$out}{$taxon_in}{$inpa};
		$co_orthologs->{$inpa}{$taxon_out}{$out} = $blastp{$out}{$taxon_in}{$inpa};
#		printf("(inserted)\t %s %s \n", $out, $inpa); printf("(inserted)\t %s %s \n", $inpa, $out);
	    } else {
		if($print_disregarded_co_orthologs) {
		    printf("(discarded_due_to_ortholog)\t %s %s \n", $out, $inpa);
		    printf("(discarded_due_to_ortholog)\t %s %s \n", $inpa, $out);
		}
	    }
	} else {
	    if($print_disregarded_co_orthologs) {
		printf("(discarded_due_to_not_found)\t %s %s \n", $out, $inpa);
		printf("(discarded_due_to_not_found)\t %s %s \n", $inpa, $out);
	    }
	}
    }
}
#! Builds the orthologs
sub build_co_orthologs {
    my ($blastp_, $inparalogs_, $orthologs_, $nss) = @_;
    my %orthologs = %$orthologs_;
# 
    my %co_orthologs;
#    printf("starts the 'build_co_orthologs(..)\n");
    foreach my $in (keys %orthologs) {
	{ #	if(defined($inparalog_limit{$in})) {
	    my $taxon_in = blast_parsing::get_taxon_type($in);
	    foreach my $taxon_out (keys %{ $orthologs{$in} }) {
		# May have several orthologs from the same taxon:
		while (my ($out, $value) = each %{ $orthologs{$in}{$taxon_out} } ) {		    		    
		    add_co_orthologs_inparalogs_ortholog($blastp_, $inparalogs_, $orthologs_, \%co_orthologs, $in, $taxon_in, $out, $taxon_out);
		    if($nss) {
			# In addition add the {(inparalogs)[taxon_a] <---> (inparalogs)[taxon_b]} relations:
			add_co_orthologs_inparalogs_inparalogs($blastp_, $inparalogs_, $orthologs_, \%co_orthologs, $in, $taxon_in, $out, $taxon_out);
		    }
		}
	    }
	} 
    }
    return %co_orthologs;
}


sub write_co_ortholog_control_file {
   my ($control_output_folder, $list_) = @_;
   blast_filtering::write_control_file($control_output_folder, $co_ortholog_name, $list_);
}

sub print_preferences {
    my ($id_2_param, $outsNotFoundInList_, $list_, $blastp_average_, $blastp_raw_) = @_;
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
    foreach my $in (keys %outsNotFoundInList) {
	my $taxon_in = blast_parsing::get_taxon_type($in);
	foreach my $taxon_out (keys %{ $outsNotFoundInList{$in} }) {
	    my $size = scalar(keys %{ $list{$in}{$taxon_out} });
	    if($size > 0) {
		printf("(**)\t Co_Orthologs for %s:\n", $in);
		while (my ($out, $value) = each %{ $outsNotFoundInList{$in}{$taxon_out} } ) {
		    blast_parsing::print_raw_data_recip($blastp_raw_, $in, $out, $taxon_in, $taxon_out, $blastp_average_);
		  }
		printf("--> \t But instead %s prefered:\n", $id_2_param);
		while (my ($out, $value) = each %{ $list{$in}{$taxon_out} } ) {
		    blast_parsing::print_raw_data_recip($blastp_raw_, $in, $out, $taxon_in, $taxon_out, $blastp_average_);
		  }
		printf("-----------------------------------------\n");
	    } else {
		printf("\t(NULL) The co_orthologs were not defined in %s, but should have been:\n", $id_2_param); 
		while (my ($out, $value) = each %{ $outsNotFoundInList{$in}{$taxon_out} } ) {
		    blast_parsing::print_raw_data_recip($blastp_raw_, $in, $out, $taxon_in, $taxon_out, $blastp_average_);
		  }
		printf("-----------------------------------------\n");
	    }
	} 
    }
}

#! Compares OrthoMcl and TurboOrtho to the control set:
sub compare_control_co_orthologs_to_other {
    my($other_id, $control_output_folder, $co_orthologs_, $other_co_orthologs_, $blastp_average_, $name_fp_other, $name_fn_other, $blastp_raw_) = @_;
    if(!(defined $blastp_raw_)) {
	printf("Not defined the blastp_raw (1); an error\n");
    }
    # The returned list holds those co_orthologs not found in orthaGogue's dataset:
    my($TP_other_1, $FN_other, %co_orthologs_fn_other) = blast_filtering::get_difference($co_orthologs_, $other_co_orthologs_);
    # Writes the lacking co_orthologs to an own file:
    blast_filtering::write_control_file($control_output_folder, $name_fn_other, \%co_orthologs_fn_other);

    # The returned list holds those co_orthologs not found in controls dataset:
    my($TP_other_2, $FP_other, %co_orthologs_fp_other) = blast_filtering::get_difference($other_co_orthologs_, $co_orthologs_);
    # Writes those errate co_orthologs to an own file:
    blast_filtering::write_control_file($control_output_folder, $name_fp_other, \%co_orthologs_fp_other);

    if($FN_other > 0) {
	printf("- The preferences are for the %d FN in %s are:\n", $FN_other, $other_id);	
	print_preferences($other_id, \%co_orthologs_fn_other, $other_co_orthologs_, $blastp_average_, $blastp_raw_);
    }
#printf("- The preferences are for the %d false hits (FP) in %s are:\n", $FP_other, $other_id);	
    if($FP_other > 0) {
	printf("- The preferences are for the %d false hits (FP) in %s are:\n", $FP_other, $other_id);	
	print_preferences("control", \%co_orthologs_fp_other, $co_orthologs_, $blastp_average_, $blastp_raw_);
    }
    if($TP_other_1 != $TP_other_2) {
	printf("- For %s, we observe that wile searcing for the FP we find %d TP or %d TP if searching for FP.\n", $other_id, $TP_other_1, $TP_other_2);
    }
    return ($TP_other_1, $FP_other, $FN_other);
}
#! Compares OrthoMcl and TurboOrtho to the control set:
sub compare_co_orthologs {
    my($control_output_folder, $co_orthologs_, $omcl_co_orthologs_, $gogue_co_orthologs_, $blastp_average_, $blastp_raw_) = @_;
     if(!(defined $blastp_raw_)) {	printf("the blastp_raw is not defined: an error in the file.\n");}
# Tests TurboOrtho:
     my($TP_gogue, $fp_gogue, $fn_gogue) = compare_control_co_orthologs_to_other("gogue", $control_output_folder, $co_orthologs_, $gogue_co_orthologs_, $blastp_average_, $name_fp_gogue, $name_fn_gogue, $blastp_raw_);
# Tests OrthoMcl:
     my($TP_omcl, $fp_omcl, $fn_omcl) = (0,0, 0);
     if((%$omcl_co_orthologs_)) {
# FIXME: Below!
 	($TP_omcl, $fp_omcl, $fn_omcl)    = compare_control_co_orthologs_to_other("omcl", $control_output_folder, $co_orthologs_, $omcl_co_orthologs_, $blastp_average_, $name_fp_omcl, $name_fn_omcl, $blastp_raw_);
     } else {printf("...\t(Does not make a control of the omcl results for the co_orthologs.)\n");}
     printf("\n#\tCo_Orthologs: In sum we see that \n-\torthaGogue has %d TP, %d FP and %d FN, while\n-\tOrthoMcl has %d TP, %d FP and %d FN.\n", $TP_gogue, $fp_gogue, $fn_gogue, $TP_omcl, $fp_omcl, $fn_omcl);
}


1;
