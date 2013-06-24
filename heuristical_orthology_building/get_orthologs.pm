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

package get_orthologs;
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
my $ortholog_name = "orthologs_control.abc";

# The output names of the differences:
my $name_fn_gogue = "orthologs_FN_gogue";
my $name_fp_gogue = "orthologs_FP_gogue";
my $name_fn_omcl = "orthologs_FN_omcl";
my $name_fp_omcl = "orthologs_FP_omcl";

# Settings changing default output:
my $print_disregarded_orthologs = 0; # If set to '1' prints the possibile orthologs discraded.
my $print_number_of_possible_orthologs = 0;

#! Builds the orthologs
sub build_orthologs {
    my ($blastp_, $inparalog_limit_, $raw_list) = @_;
    my %blastp = %$blastp_; my %inparalog_limit = %$inparalog_limit_; my %blastp_raw = %$raw_list;
    my %orthologs;
    my $cnt_poss_orth = 0;
    # Starts building a non-reciprocal ortholog set:    
    foreach my $in (keys %blastp) {
	my $taxon_in = blast_parsing::get_taxon_type($in);
	$inparalog_limit_->{$in} = 0; # If the zero-value is remained, it implies it does not have any relations in other taxa, ie, other taxa will not be able to constrain what proteins are regarded as signficant with regard to evolution.	
	foreach my $taxon_out (keys %{ $blastp{$in} }) {
	    # Only build relations for those combinations not constituting inparalogs:
	    if($taxon_in ne $taxon_out) {
		# For each taxon first find the best score:
		my $best_score = 0; 
		while (my ($out, $value) = each %{ $blastp{$in}{$taxon_out} } ) {
		    if($value > $best_score) {
			$best_score = $value;
		    }
		}
		# Then adds those having the best score as possible orthologs:
		while (my ($out, $value) = each %{ $blastp{$in}{$taxon_out} } ) {
		    if($value >= $best_score) {
			$orthologs{$in}{$taxon_out}{$out} = $value;
			$cnt_poss_orth++;
		    }
		}
		# Updates the score if it's the best:
		if($best_score > $inparalog_limit_->{$in}) {$inparalog_limit_->{$in} = $best_score;}
	    }
	}
    }
    my $cnt_orthologs = 0;
    # Converts the non-reciprocal orthologs to a reciprocal one:
    foreach my $in (keys %orthologs) {
	my $taxon_in = blast_parsing::get_taxon_type($in);
	foreach my $taxon_out (keys %{ $orthologs{$in} }) {
	    # For each taxon first find the best score:
	    my $best_score = 0; 
	    while (my ($out, $value) = each %{ $orthologs{$in}{$taxon_out} } ) {
		
		if(!(my $found = $orthologs{$out}{$taxon_in}{$in})) {
		    delete($orthologs{$in}{$taxon_out}{$out});# Deletes the ortholog.
		    } else {			
			$cnt_orthologs++;
			#! For correctness, test that the ortholog is found in the blast-file
			if(!defined($blastp_raw{$out}{$taxon_in}{$in})) {
			    printf("\t \"$out\"-->\"$in\" had no hits in the blastp-file for its reciprocal value, given taxon \"$taxon_in\", at %s:%d\n", __PACKAGE__, __LINE__);
			    exit; # TODO: remove this line!
			}
			# # FIXME: consider if below if-block is correct!
			# if(!defined($blastp_raw{$out}{$taxon_out}{$in})) {
			#     printf("...\t \"$out\"-->\"$in\" had no hits in the blastp-file for its reciprocal value, given taxon \"$taxon_in\", at %s:%d\n",__PACKAGE__, __LINE__);
			#     print("\t to test the values in the opposite direction, i.e. \"$taxon_in\" instead of \"$taxon_in\", we print these values:\n");
			#     # printf("....%s...\n", $blastp_raw{$in}{$taxon_out}{$out});
			#     # printf("....%s...\n", $blastp_raw{$out}{$taxon_in}{$in});
			#     while (my ($out_b, $value_b) = each %{ $blastp_raw{$in}{$taxon_out}{$out} } ) {
			# 	my @keys = split /\s+/, $value_b;
			# 	printf("%10s\t%10s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%-10s\n",
			# 	       $keys[0], $keys[1], $keys[2], $keys[3], $keys[4], $keys[5], $keys[6], $keys[7], $keys[8], $keys[9],
			# 	       $keys[10], $keys[11]);
			#     }
			#     while (my ($out_b, $value_b) = each %{ $blastp_raw{$out}{$taxon_in}{$in} } ) {
			# 	my @keys = split /\s+/, $value_b;
			# 	printf("%10s\t%10s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%-10s\n",
			# 	       $keys[0], $keys[1], $keys[2], $keys[3], $keys[4], $keys[5], $keys[6], $keys[7], $keys[8], $keys[9],
			# 	       $keys[10], $keys[11]);
			#     }
			#     exit; # TODO: remove this line!
			# }
		    }
	    }
	}
    }
    if($print_number_of_possible_orthologs) {    printf("#\tFound %d orthologs in the control set (of %d possible orthologs).\n", $cnt_orthologs, $cnt_poss_orth);}
#    blast_filtering::print_inparalog_limits($inparalog_limit_); # 
#    print_inparalog_limits(\%inparalog_limit); # 

    return %orthologs;
}

sub write_ortholog_control_file {
   my ($control_output_folder, $list_) = @_;
   blast_filtering::write_control_file($control_output_folder, $ortholog_name, $list_);
}

sub print_preferences {
    my ($id_2_param, $outsNotFoundInList_, $list_, $blastp_average_, $blastp_raw_) = @_;
    if(!(defined $blastp_raw_)) {
	printf("ikke definert2---en feil\n");
    }
    if(!(defined $blastp_average_)) {
	printf("'blastp_average_ not defined (1)
\n");
    }
    my %blastp_raw = %$blastp_raw_;

#    my %blastp_average = %$blastp_average_;
#    %blastp_raw = %$blastp_raw_;
    #! Internal counters:
    my $cnt_singleAlignment = 0;
    my $cnt_mult_Alignment = 0;
    my $cnt_noRecip = 0;
    # Remaps:
    my %outsNotFoundInList = %$outsNotFoundInList_; my %list = %$list_; my %blastp_average = %$blastp_average_;
    foreach my $in (keys %outsNotFoundInList) {
	my $taxon_in = blast_parsing::get_taxon_type($in);
	foreach my $taxon_out (keys %{ $outsNotFoundInList{$in} }) {
	    my $size = scalar(keys  %{ $list{$in}{$taxon_out} });
	    if($size > 0) {
		printf("\t(%d) The orthologs for %s {%s} should have been:\n", $size, $in, $taxon_out);
		while (my ($out, $value) = each %{ $outsNotFoundInList{$in}{$taxon_out} } ) {
#		    printf("#\t %s:%f\n", $out, $value);
		    blast_parsing::print_raw_data_recip($blastp_raw_, $in, $out, $taxon_in, $taxon_out, $blastp_average_);
		}
		printf("\t-> But instead %s prefered:\n", $id_2_param);
		while (my ($out, $value) = each %{ $list{$in}{$taxon_out} } ) {
#		    printf("#\t %s:%f\n", $out, $value); # printf(" (we found it having an value of %f) ", $blastp_average{$in}{$taxon_out}{$out});
		    blast_parsing::print_raw_data_recip($blastp_raw_, $in, $out, $taxon_in, $taxon_out, $blastp_average_);
		}
		printf("-----------------------------------------\n");
	    } else {
		printf("\t(NULL) The orthologs were not defined in %s, but should for %s {%s} have been:\n", $id_2_param, $in, $taxon_out);
		while (my ($out, $value) = each %{ $outsNotFoundInList{$in}{$taxon_out} } ) {
#		    printf("#\t %s:%f\n", $out, $value);#  printf(" (we found it having an value of %f) ", $blastp_average{$in}{$taxon_out}{$out});
		    blast_parsing::print_raw_data_recip($blastp_raw_, $in, $out, $taxon_in, $taxon_out, $blastp_average_);
		    #! Support evalaution of why the error occured:
		    my $size_in = scalar(keys  %{ $blastp_raw{$in}{$taxon_out}{$out}});
		    my $size_out = scalar(keys  %{ $blastp_raw{$out}{$taxon_in}{$in}});
		    #! Teast the effect of non-defined values:
		    if(($size_in == 0) && ($size_out != 0)) {
			printf("noRecip_%s \t \"$out\"-->\"$in\" had no hits in the blastp-file for its reciprocal value, given taxon \"$taxon_in\", at %s:%d\n", $id_2_param, __PACKAGE__, __LINE__);
			$cnt_noRecip++;
		    } elsif(($size_in != 0) && ($size_out == 0)) {
			printf("noRecip_%s \t \"$in\"-->\"$out\" had no hits in the blastp-file for its reciprocal value, given taxon \"$taxon_in\", at %s:%d\n", $id_2_param, __PACKAGE__, __LINE__);
			$cnt_noRecip++;
		    }		    
		    if($size_out == 1) { # then only a single alignment, implicating that the practise of choosing an HSP should have nothing to say.
			printf("singleAlignment_%s \t \"$in\"-->\"$out\" had a single match in both directions in the blastp-file, given taxon \"$taxon_in\", at %s:%d\n", $id_2_param, __PACKAGE__, __LINE__);
			$cnt_singleAlignment++;
		    } else {
			printf("mult_Alignment_%s \t \"$in\"-->\"$out\" had no hits in the blastp-file; had a single match in both directions, given taxon \"$taxon_in\", at %s:%d\n", $id_2_param, __PACKAGE__, __LINE__);
			$cnt_mult_Alignment++;
		    }
		    # if(!defined($blastp_raw{$out}{$taxon_in}{$in})) {
		    # }
		}
		printf("-----------------------------------------\n");
	    }
	} 
    }
    printf("DIFFERENCE_DESCRIPTION_%s\t The differences occured for situations for: $cnt_noRecip where reciprocal-values in blastp-file was not defined; $cnt_singleAlignment where only a single aligmentment was givne for the scores, i.e. merging of hits did not occure; $cnt_mult_Alignment for situations where the scores were summed <-- a more complex case than the first. This message was produced at %s:%d\n", $id_2_param,  __PACKAGE__, __LINE__);
}
#! Compares OrthoMcl and TurboOrtho to the control set:
sub compare_controlOrthologs_to_other {
    my($other_id, $control_output_folder, $orthologs_, $other_orthologs_, $blastp_average_, $name_fp_other, $name_fn_other, $blastp_raw_) = @_;
    if(!(defined $blastp_raw_)) {
	printf("ikke definert(1)---en feil\n");
    }
    # The returned list holds those orthologs not found in orthaGogue's dataset:
    my($TP_other_1, $FN_other, %orthologs_fn_other) = blast_filtering::get_difference($orthologs_, $other_orthologs_);
    # Writes the lacking orthologs to an own file:
    blast_filtering::write_control_file($control_output_folder, $name_fn_other, \%orthologs_fn_other);

    # The returned list holds those orthologs not found in controls dataset:
    my($TP_other_2, $FP_other, %orthologs_fp_other) = blast_filtering::get_difference($other_orthologs_, $orthologs_);
    # Writes those errate orthologs to an own file:
    blast_filtering::write_control_file($control_output_folder, $name_fp_other, \%orthologs_fp_other);

    if($FN_other > 0) {
	printf("- The preferences are for the %d FN in %s are:\n", $FN_other, $other_id);	
	my $control = "FN__" . $other_id;
	print_preferences($control, \%orthologs_fn_other, $other_orthologs_, $blastp_average_, $blastp_raw_);
    }
#printf("- The preferences are for the %d false hits (FP) in %s are:\n", $FP_other, $other_id);	
    if($FP_other > 0) {
	printf("- The preferences are for the %d false hits (FP) in %s are:\n", $FP_other, $other_id);	
#	print_preferences("control", \%orthologs_fp_other, $other_orthologs_, $blastp_average_, $blastp_raw_); 
	my $control = "FP__" . $other_id;
	print_preferences($control, \%orthologs_fp_other, $orthologs_, $blastp_average_, $blastp_raw_); 
    }
    if($TP_other_1 != $TP_other_2) {
	printf("- For %s, we observe that wile searcing for the FP we find %d TP or %d TP if searching for FP.\n", $other_id, $TP_other_1, $TP_other_2);
    }
    return ($TP_other_1, $FP_other, $FN_other);
}
#! Compares OrthoMcl and TurboOrtho to the control set:
sub compare_orthologs {
    my($control_output_folder, $orthologs_, $omcl_orthologs_, $gogue_orthologs_, $blastp_average_, $blastp_raw_) = @_;
    if(!(defined $blastp_raw_)) {	printf("the blastp_raw is not defined: an error in the file.\n");}
# Tests TurboOrtho:
    my($TP_gogue, $fp_gogue, $fn_gogue) = compare_controlOrthologs_to_other("orthAgogue", $control_output_folder, $orthologs_, $gogue_orthologs_, $blastp_average_, $name_fp_gogue, $name_fn_gogue, $blastp_raw_);
# Tests OrthoMcl:
    my($TP_omcl, $fp_omcl, $fn_omcl) = (0,0, 0);
    if((%$omcl_orthologs_)) {
	($TP_omcl, $fp_omcl, $fn_omcl)    = compare_controlOrthologs_to_other("omcl", $control_output_folder, $orthologs_, $omcl_orthologs_, $blastp_average_, $name_fp_omcl, $name_fn_omcl, $blastp_raw_);
    } else {printf("...\t(Does not make a control of the omcl results for the orthologs.)\n");}
#    my($TP_omcl, $fp_omcl, $fn_omcl)    = compare_controlOrthologs_to_other("omcl", $control_output_folder, $orthologs_, $omcl_orthologs_, $blastp_average_, $name_fp_omcl, $name_fn_omcl, $blastp_raw_);
    printf("\n#\tOrthologs: In sum we see that \n-\torthAgogue has %d TP, %d FP and %d FN, while\n-\tOrthoMcl has %d TP, %d FP and %d FN.\n", $TP_gogue, $fp_gogue, $fn_gogue, $TP_omcl, $fp_omcl, $fn_omcl);

}


1;
