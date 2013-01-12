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

package blast_parsing;
use warnings;
use strict;
use POSIX;
#use types_and_enums;
#use types_and_enums qw($format_blastp_evalue $format_blastp_last $format_abc );
# Defines the booleans in the code, making it easier to shift this code to c++ (as always is or path of walk).
use constant false => 0;
use constant true  => 1;
#package blast_parsing;
# Defines the booleans in the code, making it easier to shift this code to c++ (as always is or path of walk).
use constant false => 0;
use constant true  => 1;
# # Defines the file format to parse:
use constant format_blastp_evalue   => 'blastp_c_11';
use constant format_blastp_last     => 'blastp_c_12';
use constant format_abc             => 'abc';

# Global variables for identifying the taxa:
my $taxon_index = 0;
my $taxon_protein_split = '\|'; # Note the mandatory 'dash', as the pipe is a special symbol!
my $NOT_append_data = 0;

# Local settings updated during the process:
#my $maximum_distance = 0; # The maxium distance;
# FIXME: Include the above, and remove the below!
my $maximum_distance = 10; # The maxium distance;

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
    my @taxon_1 = get_taxon_type($label_tot[0]);
    my @taxon_2 = get_taxon_type($label_tot[1]);
    return ($taxon_1[$taxon_index], $taxon_2[$taxon_index]);
}

#! Sets the maximum distance;
sub update_maximum_distance {
    my ($distance) = @_;
    if(defined $distance) {
	if(isdigit($distance) && $distance > $maximum_distance) {
	    $maximum_distance = $distance;
	}
    }
}

#! Returns the distance
sub get_distance {
    my ($format, @keys) = @_;
    if((defined $format)) { # && (scalar(@keys)>0) {
	my $distance = 0;
	if($format eq format_abc) {
	    $distance = $keys[2]; 
	} else {
	    if($format eq format_blastp_evalue) {
		if($keys[10] == 0) {
		    $distance = $maximum_distance +1;
#		printf("\tkeys[10] = %e && distance=%f+1\n", $keys[10], $maximum_distance);
		} else {
		    $distance = -1*log10($keys[10]); 	    
#		printf("\tkeys[10] = %e && distance=%f\n", $keys[10], $distance);
		}
	    } else {
		$distance = $keys[11]; 
	    }
	}
	return $distance; 
    } else {
	return 0;
    }
}
#! Return true if distance is to be appended to the repset value.
sub get_appended_distance {
    my ($format, $new_distance, $old_distance) = @_;
    my $correct_distance = $new_distance;
    if($format ne format_abc) { # FIXME!
	if(defined $old_distance && ($old_distance > 0)) { # a value exists:
	    if(0==$NOT_append_data) { #NOT(NOT) gives "append data".
		$correct_distance= ($old_distance+$new_distance);
	    } else {
		if($old_distance > $new_distance) { $correct_distance= $old_distance;}
	    } 
	} 
    }
#    printf("#\treturns %d, given old_distance(%f) and new_distance(%f)\n", $correct_distance, $old_distance, $new_distance);
    return $correct_distance;
} 

##return $new_distance;}
# sub shall_append_distance {
#     my ($format, $current, $new_distance) = @_;
#     if($format ne format_abc) { # FIXME!
# 	if(defined $current) {
# 	    if(!$NOT_append_data) { #NOT(NOT) gives "append data".
# 		return true;
# 	    } else {return false;}
# 	}
#     } 
#     return false;
# }
#! Sets the maximum value and produces the list of taxa:
sub get_blastp_info {
    my ($file, $format) = @_;
    if(-e $file) {
	open(my $file1, "<", $file) or die $!;
	my %taxa;
	while (<$file1>) {	
	    my @keys = split /\s+/, $_;
	    my ($taxon1, $taxon2) = get_taxa(@keys);
# FIXME: remove the printf below!
#	    printf("Has format=%s and keys=%s, in blast_parsing.pm\n", $format, @keys);
	    my $distance = get_distance($format, @keys);
	    update_maximum_distance($distance);
	    if(defined $taxa{$taxon1}{$taxon2}) {
		$taxa{$taxon1}{$taxon2}++;
	    } else {
		$taxa{$taxon1}{$taxon2} = 1;
	    }	
	}
	close($file1) or die $!;
	return %taxa;
    } else {
	printf("!!\tNot able to fined blastp file: %s\n", $file);
	return NULL;
    }
}
# sub trim($)
# {
#     my $string = shift;
#     $string =~ s/\s+$//; # removes trailing white spaces
# #    $string =~ s/\s+$//; # Removes new-lines
# #    $string =~ s/^\s+//;
# #    $string =~ s/\s+$//;
#     return $string;
# }
sub print_blastp_data {
    my(@keys) = @_;
    my $size = scalar(@keys);
    my $cnt = 1;
    {
	my $i = 0;
	while($i < $size-2) {
	    print "\$keys[$i], ";
#	    printf("%s\t", $keys[$i]);
	    $i++;
	}
#	printf("%-10s", $keys[$i]);
#	printf("%s ", $keys[$i+1]);
    }
}
#! Prints the raw data found in the lbastp data for the given proteins.
sub print_raw_data_recip {
    my ($raw_list, $in, $out, $taxon_in, $taxon_out, $blastp_average_, $limit_in, $limit_out) = @_;
    my %blastp_raw = %$raw_list;
    my %blastp_average = %$blastp_average_;
    if(!(defined $blastp_average_)) {
	printf("'blastp_average_ not defined (2)\n");
    }
    my $size_in = scalar(keys  %{ $blastp_raw{$in}{$taxon_out}{$out}});
    my $size_out = scalar(keys  %{ $blastp_raw{$out}{$taxon_in}{$in}});
    if($size_in > 0 || $size_out > 0) { 
	# Note: In order to show the algortihms procedure, the 'transformed'- and averaged values are included on the each line, enclosed in paranthesis at the lines end.
	while (my ($out_b, $value_b) = each %{ $blastp_raw{$in}{$taxon_out}{$out} } ) {
	    my @keys = split /\s+/, $value_b;
	    chomp($value_b); # Removes the newline char.
	    my $value_avg = 0; my $avg_poss = $blastp_average{$in}{$taxon_out}{$out};
	    if(defined($avg_poss)) { $value_avg = $avg_poss;} #$blastp_average{$in}{$taxon_out}{$out_b};}
	    my $distance = get_distance(format_blastp_evalue, @keys);
	    if(defined($blastp_average{$in}{$taxon_out}{$out_b})) {$value_avg=$blastp_average{$in}{$taxon_out}{$out_b};}
	    printf("%10s\t%10s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%-10s  (%f  %f",
		   $keys[0], $keys[1], $keys[2], $keys[3], $keys[4], $keys[5], $keys[6], $keys[7], $keys[8], $keys[9],
		   $keys[10], $keys[11], $distance, $value_avg
		   );
	    if(defined $limit_in && defined $limit_out) {printf(" limit_min(%f, %f)", $limit_out, $limit_in);}

	    printf(")\n");
	}
#	printf("                                                                                                                                     ...................\n");
	while (my ($out_b, $value_b) = each %{ $blastp_raw{$out}{$taxon_in}{$in} } ) {
	    my @keys = split /\s+/, $value_b;
	    chomp($value_b);
	    my $value_avg = 0; my $avg_poss = $blastp_average{$in}{$taxon_out}{$out};
	    if(defined($avg_poss)) { $value_avg = $avg_poss;} #$blastp_average{$in}{$taxon_out}{$out_b};}
	    my $distance = get_distance(format_blastp_evalue, @keys); 
	    printf("%10s\t%10s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%-10s  (%f  %f",
		   $keys[0], $keys[1], $keys[2], $keys[3], $keys[4], $keys[5], $keys[6], $keys[7], $keys[8], $keys[9],
		   $keys[10], $keys[11], $distance, $value_avg
		   );
	    if(defined $limit_in && defined $limit_out) {printf(" limit_min(%f, %f)", $limit_out, $limit_in);}
	    printf(")\n");
	}
    }
}
#! Returns the averaged list:
sub make_avarage {
    my ($list_, $blastp_raw_) = @_;
    my %list = %$list_; my %blastp_raw = %$blastp_raw_;
    my %average;
    foreach my $in (keys %list) {
	foreach my $taxon_out (keys %{ $list{$in} }) {
	    while (my ($out, $value) = each %{ $list{$in}{$taxon_out} } ) {
		my $taxon_in = blast_parsing::get_taxon_type($in);
		if (!(my $found = $list{$out}{$taxon_in}{$in})) {
	#	    delete($list{$in}{$taxon_out}{$out}); # Removes it.
		    # The above line changed the result of the co-orthologs: If not fully understood, comment out the above, and comment in the below line:
		    $average{$out}{$taxon_in}{$in} = $average{$in}{$taxon_out}{$out} = $value/2;
#$average{$out}{$taxon_in}{$in} = $value/2;
		} else { # Sets the average:
		    my $average_value = ($list{$out}{$taxon_in}{$in} + $list{$in}{$taxon_out}{$out})/2;;
		    $average{$out}{$taxon_in}{$in} = $average{$in}{$taxon_out}{$out} = $average_value;
		}
	    }
	}
    }    
    return %average;
}


#! Gets the hash itself (from the file):
sub get_hash {
    my ($file, $format, $NOT_append_data_) = @_;
    $NOT_append_data = $NOT_append_data_;
#    printf("tries open file-name '%s'\n", $file);
    if(-e $file) {
	open(my $file1, "<", $file); # or die $!;
	my %file1keys;
	my $cnt_new_keys_inserted = 0;
	while (<$file1>) {	
	    my @keys = split /\s+/, $_;
#	my $taxon2 = get_taxon(@keys);
#	my $str = $_;	chomp($str);	printf("%s\n", $str);
	    my ($taxon1, $taxon2) = get_taxa(@keys);
	    my $new_distance = get_distance($format, @keys);
	    if((defined $keys[0]) && (defined $taxon2) && (defined $keys[1])) {
#	    if(defined $file1keys{$keys[0]}{$taxon2}{$keys[1]}) {
		my $old_distance = $file1keys{$keys[0]}{$taxon2}{$keys[1]};
#	if(defined($old_distance)) {	    printf("Tests if %f > %f\n", $old_distance, $new_distance);	}
		
		my $new_distance_append = get_appended_distance($format, $new_distance, $old_distance);
		if(defined $new_distance) {
		    $file1keys{$keys[0]}{$taxon2}{$keys[1]} = $new_distance_append;
		}
	    }
#	printf("value set to %f (new_distance was %d)\n", 	$file1keys{$keys[0]}{$taxon2}{$keys[1]}, $new_distance);
# 	if(!shall_append_distance($format, $distance)) {
# 	    $cnt_new_keys_inserted++;	    
# 	}
#	$file1keys{$keys[0]}{$taxon2}{$keys[1]} = $distance;
# 	if(!shall_append_distance($format, $distance)) {
# 	    $file1keys{$keys[0]}{$taxon2}{$keys[1]} = $distance;
# 	    $cnt_new_keys_inserted++;	    
# 	} else {
# 	    $file1keys{$keys[0]}{$taxon2}{$keys[1]} += $distance;
# 	}
#	printf("---------------------\n");
	}
	close($file1) or die $!;
    return %file1keys;
    }
#    printf("From the file %s found in total %d unique pairs\n", $file, $cnt_new_keys_inserted);
    return ();
}

#! Puts the raw blastp file into a hash for easy access later on:
sub get_raw_blastp_file {
    my ($file, $format) = @_;
    open(my $file1, "<", $file) or die $!;
    my %file1keys;
    while (<$file1>) {	
 	my @keys = split /\s+/, $_;
#	my $taxon2 = get_taxon(@keys);
	my ($taxon1, $taxon2) = get_taxa(@keys);	
	my $distance = get_distance($format, @keys);
	my $size = scalar(keys  %{ $file1keys{$keys[0]}{$taxon2}{$keys[1]} });
	# Inserts at the end of it:
	$file1keys{$keys[0]}{$taxon2}{$keys[1]}{$size} = $_; 
	$size = scalar(keys  %{ $file1keys{$keys[0]}{$taxon2}{$keys[1]} });
#	printf("found protein: %s %s (with innter taxon %s and outer taxon %s and size %d)\n", $keys[0], $keys[1], $taxon1, $taxon2, $size);
    }
    close($file1) or die $!;
    return %file1keys;
}


sub print_taxa {
    my ($taxa_) = @_;
    my %taxa = %$taxa_;
    printf("\nStudies the set of taxon-taxon pairs with notation <outer taxon>:<number of protein pairs>:\n");
    foreach my $in (keys %taxa) {
	printf("- Looks at taxon-clusters for '%s': ", $in);
	while (my ($out, $value) = each %{ $taxa{$in} } ) {
	    printf("%s:%d ", $out, $value);
	}
	printf("\n");
    }
    printf("------------------------------------------------------\n");
}
1;
