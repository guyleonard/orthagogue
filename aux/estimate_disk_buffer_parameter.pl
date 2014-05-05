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

# # Below a trick forgotten for a c++ programmer (oekseth)
# use File::Basename;
# use lib basename ($0);
# use FindBin;
# use lib "$FindBin::Bin";
# use File::Copy;

# Change [below] in order to increase the memory consumption.
my %settings;
$settings{"update_taxa_count"} = 0;


# Defines the booleans in the code, making it easier to shift this code to c++ (as always is or path of walk).
use constant false => 0;
use constant true  => 1;
# Defines the booleans in the code, making it easier to shift this code to c++ (as always is or path of walk).
use constant false => 0;
use constant true  => 1;
# Defines the file format to parse:
use constant format_blastp_evalue   => 'blastp_c_11';
use constant format_blastp_last     => 'blastp_c_12';
use constant format_abc             => 'abc';

# Global variables for identifying the taxa:
my $taxon_index = 0;
my $taxon_protein_split = '\|'; # Note the mandatory 'dash', as the pipe is a special symbol!

my $biggest_identified_disk_buffer_size_taxon_taxon = 0;
my $biggest_identified_disk_buffer_size_taxon_taxon_taxaPair = "";
my $biggest_identified_disk_buffer_size_taxon_taxon_protein = 0;
my $biggest_identified_disk_buffer_size_taxon_taxon_protein_id = "";
my $biggest_identified_disk_buffer_size_taxon_taxon_protein_taxaPair = "";



#! Returns the taxon type of the given argument.
sub get_taxon_type {
    my($label) = @_;
    if(defined $label) {
	my @taxon = split($taxon_protein_split, $label);
	return $taxon[$taxon_index];
    } else {
	return 0;
    }
}
#! Returns the taxa using tha above global variables:
sub get_taxa {
    my (@label_tot) = @_;
    my $taxon_1 = get_taxon_type($label_tot[0]);
    my $taxon_2 = get_taxon_type($label_tot[1]);
    return ($taxon_1, $taxon_2);
}


sub update_taxa_count {
    my ($taxa, $keys) = @_;
    if(scalar(@_) != 2) {warn("Wrong number of arguments provided to function.\n");}
    my ($taxon1, $taxon2) = get_taxa(@{$keys}); 
    if(!defined($taxon1)) {
	warn("taxon1 not defined, given taxon-protein-seperator=\"$taxon_protein_split\" and taxon-index=\"$taxon_index\" at [%s]:%d\n", __PACKAGE__, __LINE__);
	($taxon1, $taxon2) = (undef, undef);
	exit; # todo: consider removing this line.
    } elsif(!defined($taxon2)) {
	warn("taxon2 not defined, given taxon-protein-seperator=\"$taxon_protein_split\" and taxon-index=\"$taxon_index\" at [%s]:%d\n", __PACKAGE__, __LINE__);
	($taxon1, $taxon2) = (undef, undef);
	exit; # todo: consider removing this line.
    } else {
	if(defined $taxa->{$taxon1}{$taxon2}) {
	    $taxa->{$taxon1}{$taxon2}++;
	} else {
	    $taxa->{$taxon1}{$taxon2} = 1;
	}	
    }
    return ($taxon1, $taxon2);
}


#! Get the number of pairs, and number of chars, for [taxa][taxa][protein].
sub update_inner_protein_count {
    my ($protein1, $list_taxa_protein_cnt, $list_taxa_protein_charUsage,  $cnt_chars_in_line, $taxon1, $taxon2) = @_;
    if(scalar(@_) != 6) {warn("Wrong number of arguments provided to function.\n");}

    if(defined($list_taxa_protein_cnt->{$taxon1}{$taxon2}{$protein1})) {
	$list_taxa_protein_cnt->{$taxon1}{$taxon2}{$protein1}++;
    } else {
	$list_taxa_protein_cnt->{$taxon1}{$taxon2}{$protein1} = 1;
    }	
    if(defined($list_taxa_protein_charUsage->{$taxon1}{$taxon2}{$protein1})) {
	$list_taxa_protein_charUsage->{$taxon1}{$taxon2}{$protein1} += $cnt_chars_in_line;
    } else {
	$list_taxa_protein_charUsage->{$taxon1}{$taxon2}{$protein1} = $cnt_chars_in_line;
    }
    #! Then identify the biggest disk-buffer-size required (to encapsulate a complete taxa-taxa-protein-pair):
    if($list_taxa_protein_charUsage->{$taxon1}{$taxon2}{$protein1} > $biggest_identified_disk_buffer_size_taxon_taxon_protein) {
	$biggest_identified_disk_buffer_size_taxon_taxon_protein = $list_taxa_protein_charUsage->{$taxon1}{$taxon2}{$protein1};
	$biggest_identified_disk_buffer_size_taxon_taxon_protein_id = $protein1;
	$biggest_identified_disk_buffer_size_taxon_taxon_protein_taxaPair = $taxon1 . "<-->" . $taxon2;
    }
}

#! @return the number of chars is the given string.
sub get_number_of_chars {
    my ($string) = @_;
    if(scalar(@_) != 1) {warn("Wrong number of arguments provided to function.\n");}
    my @chars = split("", $string);
    return scalar(@chars);
}



#! Sets the maximum value and produces the list of taxa:
sub get_blastp_info {
    my ($file, $format, $list_taxa_charUsage, $list_taxa_protein_cnt, $list_taxa_protein_charUsage) = @_;
    if(scalar(@_) != 5) {warn("Wrong number of arguments provided to function.\n");}
    if(-e $file) {
	open(my $file1, "<", $file) or die $!;
	my %taxa;
	while (<$file1>) {	
	    my $string = $_;
	    my @keys = split /\s+/, $string;
	    my ($taxon1, $taxon2) = update_taxa_count(\%taxa, \@keys);
	    if(defined($taxon1)) {
		my $char_count = get_number_of_chars($string);
		if(defined($list_taxa_charUsage->{$taxon1}{$taxon2})) {
		    $list_taxa_charUsage->{$taxon1}{$taxon2} += $char_count;
		} else {
		    $list_taxa_charUsage->{$taxon1}{$taxon2} = $char_count;
		}
		#! Then identify the biggest disk-buffer-size required (to encapsulate a complete taxa-taxa-pair):
		if($list_taxa_charUsage->{$taxon1}{$taxon2} > $biggest_identified_disk_buffer_size_taxon_taxon) {
		    $biggest_identified_disk_buffer_size_taxon_taxon = $list_taxa_charUsage->{$taxon1}{$taxon2};
		    $biggest_identified_disk_buffer_size_taxon_taxon_taxaPair = $taxon1 . "<-->" . $taxon2;
		}
		if((defined($settings{"update_protein_count"}) && $settings{"update_protein_count"} == 1)) {
		    #! Get the number of pairs, and number of chars, for [taxa][taxa][protein].
		    update_inner_protein_count($keys[0], $list_taxa_protein_cnt, $list_taxa_protein_charUsage, $char_count, $taxon1, $taxon2);
		}
	    }

	}
	close($file1) or die $!;
	return %taxa;
    } else {
	printf("!!\tNot able to fined blastp file: %s\n", $file);
	return NULL;
    }
}

#! Outprint the cluster size of the [taxon][taxon] pairs.
sub print_taxa {
   my ($taxa, $list_taxa_charUsage) = @_;
   if(scalar(@_) != 2) {warn("Wrong number of arguments provided to function.\n");}
   my $file_name = "taxa_facts.txt";
   open(FILE, ">", $file_name) or die("Unable to open the file $file_name\n");

    printf("\n The [taxon][taxon] cluster size (with notation of \"<outer taxon>:(<number of protein pairs>,disk-buffer-size)\"):\n");
    printf(FILE "\n The [taxon][taxon] cluster size (with notation of \"<outer taxon>:(<number of protein pairs>,disk-buffer-size)\"):\n");

    foreach my $in (keys %{$taxa}) {
	printf("- Looks at taxon-clusters for '%s': ", $in);
	printf(FILE "- Looks at taxon-clusters for '%s': ", $in);
	while (my ($out, $value) = each %{ $taxa->{$in} } ) {
	    my $disk_buffer_size = $list_taxa_charUsage->{$in}{$out};
	    my $disk_buffer_size_MB = $disk_buffer_size / (1024*1024);
#	    my $disk_buffer_size_string = $disk_buffer_size . "B";
	    my $disk_buffer_size_string_MB = $disk_buffer_size_MB . "MB";
	    printf("%s:(%d,$disk_buffer_size_string_MB) ", $out, $value);
	    printf(FILE "%s:(%d,$disk_buffer_size_string_MB) ", $out, $value);
	}
	printf("\n");
	printf(FILE "\n");
    }
    printf("------------------------------------------------------\n");
    printf(FILE "------------------------------------------------------\n");
    printf(FILE "-\t In brief, w.r.t. the disk_buffer_size parameter, we observed the biggest taxon-taxon-pair consists of %u=%.3E chars (for taxon-pair \"$biggest_identified_disk_buffer_size_taxon_taxon_taxaPair\"), while\n", $biggest_identified_disk_buffer_size_taxon_taxon, $biggest_identified_disk_buffer_size_taxon_taxon);
    #! Close the file:
    close(FILE);
}

#! Outprint the cluster size of the [taxon][taxon] pairs.
sub print_inner_protein_facts {
    if((defined($settings{"update_protein_count"}) && $settings{"update_protein_count"} == 1)) {
	my ($list_taxa_protein_cnt, $list_taxa_protein_charUsage) = @_;
	if(scalar(@_) != 2) {warn("Wrong number of arguments provided to function.\n");}
	
	my $file_name = "protein_facts.txt";
	open(FILE, ">", $file_name) or die("Unable to open the file $file_name\n");
	
	printf("\n The [taxon][taxon] cluster size (with notation of \"<outer taxon>:<number of protein pairs>\"):\n");
	printf(FILE "\n The [taxon][taxon] cluster size (with notation of \"<outer taxon>:<number of protein pairs>\"):\n");
	
	foreach my $in (keys %{$list_taxa_protein_cnt}) {
	    while (my ($out, $value) = each %{ $list_taxa_protein_cnt->{$in} } ) {
		printf("- Looks at protein-clusters for taxon-pair '[$in][$out]':\n");
		printf(FILE "- Looks at protein-clusters for taxon-pair '[$in][$out]':\n");
		while (my ($protein1, $cnt_outer_proteins) = each %{ $list_taxa_protein_cnt->{$in}{$out} } ) {
		    my $disk_buffer_size = $list_taxa_protein_charUsage->{$in}{$out}{$protein1};
		    my $disk_buffer_size_MB = $disk_buffer_size / (1024*1024);
		    my $disk_buffer_size_string = $disk_buffer_size . "B";
		my $disk_buffer_size_string_MB = $disk_buffer_size_MB . "MB";
		    printf("\t protein \"$protein1\" has \"$cnt_outer_proteins\" outer proteins, and requires a disk_buffer_size=\"$disk_buffer_size_string\"=\"$disk_buffer_size_string_MB\".\n");
		    printf(FILE "\t protein \"$protein1\" has \"$cnt_outer_proteins\" outer proteins, and requires a disk_buffer_size=\"$disk_buffer_size_string\"=\"$disk_buffer_size_string_MB\".\n");
		}
	    }
	    printf("\n");
	    printf(FILE "\n");
	}
    } else {
	printf("\t The protein-count-option was de-activated in order to reduce memory consumption\n");
    }

   printf(FILE "-\t in brief, w.r.t. the disk_buffer_size parameter, we observed for the biggest taxon-taxon-protein-pair, consists of %u=%.3E chars (for protein \"$biggest_identified_disk_buffer_size_taxon_taxon_protein_id\" and taxon-pair \"$biggest_identified_disk_buffer_size_taxon_taxon_protein_taxaPair\")\n", $biggest_identified_disk_buffer_size_taxon_taxon_protein, $biggest_identified_disk_buffer_size_taxon_taxon_protein);

    printf("------------------------------------------------------\n");
    printf(FILE "------------------------------------------------------\n");

    #! Close the file:
    close(FILE);
}

sub control_files {
    my($blastp_file) = @_;
    if(scalar(@_) != 1) {warn("Wrong number of arguments provided to function.\n");}
#
# The parsing of the blastp file used for both omcl and to:
    printf("\t Starts analyzing the \"$blastp_file\", with taxon-protein-seperator=\"$taxon_protein_split\" and taxon-index =\"$taxon_index\" at index blastP file.\n");
    my %list_taxa_charUsage;
    my %list_taxa_protein_charUsage;
    my %list_taxa_protein_cnt;
    my %taxa = get_blastp_info($blastp_file, format_blastp_evalue, \%list_taxa_charUsage, \%list_taxa_protein_cnt, \%list_taxa_protein_charUsage);
    
    #! Always remember the maximum number of [taxa][taxa] and [taxa][taxa][protein], and write this result to both the file, and the stdout 'channel'.
    printf("# In brief, w.r.t. the disk_buffer_size parameter, we observed:\n");
    printf("-\t the biggest taxon-taxon-pair consists of %u=%.3E chars (for taxon-pair \"$biggest_identified_disk_buffer_size_taxon_taxon_taxaPair\"), while\n", $biggest_identified_disk_buffer_size_taxon_taxon, $biggest_identified_disk_buffer_size_taxon_taxon);
    printf("-\t the biggest taxon-taxon-protein-pair consists of %u=%.3E chars (for protein \"$biggest_identified_disk_buffer_size_taxon_taxon_protein_id\" and taxon-pair \"$biggest_identified_disk_buffer_size_taxon_taxon_protein_taxaPair\")\n", $biggest_identified_disk_buffer_size_taxon_taxon_protein, $biggest_identified_disk_buffer_size_taxon_taxon_protein);
    
    print_taxa(\%taxa, \%list_taxa_charUsage);
    print_inner_protein_facts(\%list_taxa_protein_cnt, \%list_taxa_protein_charUsage);
}


my $numArgs = $#ARGV + 1;
if(($numArgs == 3) || ($numArgs == 4)) {
    #! Update the sepcification:
    my $get_protein_facts = 0;
    $taxon_protein_split = $ARGV[1];     $taxon_index = $ARGV[2]; $get_protein_facts = $ARGV[3];
    if(defined($get_protein_facts) && ($get_protein_facts == 1) ) {
	$settings{"update_protein_count"} = 1;
    } else {
	$settings{"update_protein_count"} = 0;
    }
    #! Then make 'the call':
    control_files($ARGV[0]);
} else {
    printf("----------------------------------------------------------------------\n");
    printf("Extract BlastP facts, which (among others) is of interest when configuring orthAgogue for large blastP files.\n");
    printf("-\tUsage: %s <blastp_file> <seperator> <taxon-index> <bool: get-expressive-knowledge-of-proteins> \n", $0);
    printf("-->\tThis message was seen because %d arguments were used. In order to run the software, provide either 1 argument (i.e. the blastp-file-name).\n", $numArgs);
    printf("\nThe software was developed by O.K. Ekseth. Questions to be forwarded to [oekseth\@gmail.com].\n");
    printf("----------------------------------------------------------------------\n");
}
