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

package blast_filtering;
use warnings;
use strict;
use POSIX;
#use types_and_enums qw($format_blastp_evalue $format_blastp_last $format_abc );
use blast_parsing;
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

# Names of global variables in use:
#my %blastp_average;
# Names of the files generated:
my $co_orthologs_name = "co_orthologs_control.abc";
my $name_raw_averaged_data = "averaged_data.abc";



#! Consistency script ensuring that the building of the folder is done.
sub make_folder_if_not_set {
    my ($control_output_folder) = @_;
    my $retval = open(my $file1, "<", $control_output_folder); # or die $!;
    if(defined $retval) {
	close($file1) or die $!;
	return 1;
    } else {
	mkdir $control_output_folder, 0755;
	return 0;
    }
}

sub write_control_file {
   my ($control_output_folder, $file_name, $list_) = @_;
   my %list = %$list_;
   make_folder_if_not_set($control_output_folder);
   open(my $file_out, ">", $control_output_folder . $file_name) or die $!;
    foreach my $in (keys %list) {
	foreach my $taxon_out (keys %{ $list{$in} }) {
	    while (my ($out, $value) = each %{ $list{$in}{$taxon_out} } ) {
		printf($file_out "%s\t%s\t%f\n", $in, $out, $value);
	    }
	}
    }

   close($file_out) or die $!;
}

#! Comparing first with second hash:
sub compare_verbouse {
    my ($id, $file1keys_, $file2keys_) = @_;
    my %file1keys = %$file1keys_; my %file2keys = %$file2keys_;  # Dereferencing them for easier access.
    my $cnt = 0;
    foreach my $in (keys %file1keys) {
#	my $taxon_in = get_taxon_type($in);
	foreach my $taxon_out (keys %{ $file1keys{$in} }) {
#	while (my ($out, $value) = each %{ $file1keys{$in} } ) {
	    while (my ($out, $value) = each %{ $file1keys{$in}{$taxon_out} } ) {
		if (!(my $found = $file2keys{$in}{$taxon_out}{$out})) {
		    printf("%s\tNot found '%s %s'\n", 		       $id, $in, $out);
		    $cnt++;
		} #else {printf("%s %s\n", $in,  $out);}
	    }
	}
   }
#    printf("------------------------------------------------------\n");
    return $cnt;
}

#! Comparing first with second hash:
sub get_difference {
    my ($file1keys_, $file2keys_) = @_;
    # Remaps:
    my %file1keys = %$file1keys_; my %file2keys = %$file2keys_;  # Dereferencing them for easier access.
    my $cnt_equal = 0; # The number of htis the second argument corresponds to the first;
    my $cnt_difference = 0; # The number of htis the second argument differs from the first;
    my %difference = ();
    foreach my $in (keys %file1keys) {
	foreach my $taxon_out (keys %{ $file1keys{$in} }) {
	    while (my ($out, $value) = each %{ $file1keys{$in}{$taxon_out} } ) {
		if (!(my $found = $file2keys{$in}{$taxon_out}{$out})) {
		    $cnt_difference++;
		    $difference{$in}{$taxon_out}{$out} = $value;
		} else {
		    $cnt_equal++;		    
		}
	    }
	}
    }
    return ($cnt_equal, $cnt_difference, %difference);
}



#! Writes the averaged data used as basis for the operation.
sub write_averaged_data_file {
    my($control_output_folder, $list_) = @_;
    write_control_file($control_output_folder, $name_raw_averaged_data, $list_);
}

sub print_inparalog_limits {
    my ($inparalog_limit_) = @_;
    my %inparalog_limit = %$inparalog_limit_;
    printf("The inparalog-limits are:\n");
    foreach my $in (keys %inparalog_limit) {
	printf("%s\t%f\n", $in, $inparalog_limit{$in});	
    }
}

1; # The return value
