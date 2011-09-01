#
#===============================================================================
#
#         FILE:  Nexus_partitioning.pm
#
#  DESCRIPTION:  This module is used to detect if a partition scheme is indicated 
#                in a nexus alignment file. If exists, it can create different
#                aligments files in the desired format using bioperl module 
#                AlignIO. 
#
#        FILES:  ---
#         BUGS:  ---
#        NOTES:  Bioperl modules are needed.
#       AUTHOR:   (Francesc Lopez-Giraldez), <francesc.lopez@gmail.com>
#      COMPANY:  Yale University
#      VERSION:  0.1
#      CREATED:  02/16/2009 04:06:25 PM EST
#     REVISION:  ---
#===============================================================================

package Nexus_partitioning;

use strict;
use warnings;
use Carp;
use Bio::AlignIO;

#Global variables
my $path_to_user_folder;
my $nexus_file;
my $align_object;
my $output_format;
my $file_content;
my $align_length;
my $begin_set;

#######################################
# Partitioning the Nexus alignment in 
# slices.if the partition statement 
# is found in the nexus file.
#######################################

sub get_partitions {
    ($nexus_file)          = shift @_;
    ($align_object)        = shift @_;
    ($output_format)       = shift @_;
    $align_length          = $align_object->length();

    #Getting the file content in a string
    open  my $NEXUSFILE, '<', "$nexus_file"
        or croak "$0 : failed to open  input file '$nexus_file' : $!\n";
    
    #slurp file
    $file_content = do { local $/; <$NEXUSFILE> };

    close  $NEXUSFILE
        or carp "$0 : failed to close input file '$nexus_file' : $!\n";

    # checking if partitions are set
    # if they are set, the format is checking step by step
    # with a series of checkpoint for errors.
    my $partition_name = are_partition_sets();
    return 0 if (!$partition_name);
    
    ($begin_set) = $file_content =~ /^\s+(begin\s+sets.+^\s+end)/ixms;
    
    my ($partition_number, @partition_names) = get_partition_names($partition_name);

    if (!$partition_number || !@partition_names) {
        return 0, 1, "No partition name list was found for [$partition_name].";
    }
    elsif ($partition_number != scalar @partition_names) {
        return 0, 1, "# of partitions [$partition_number] different from # of partition names["
                         . scalar @partition_names . '].';
    }

    my ($ok, $msg, $partition_coordinates_for_ref) = get_partition_coordinates(@partition_names);
    return 0, 1, "$msg" if ($ok == 0);

    my %partition_coordinates_for = %{$partition_coordinates_for_ref};

    my $num_charset_found = scalar keys %partition_coordinates_for;
    if ($partition_number != $num_charset_found) {
        return 0, 1, "# of partitions [$partition_number] different from"
                            . " # of charset found [$num_charset_found].";
    }

    #Here we actually create the partitions, once all checkpoint are ok
    #each partition will have its own file
    create_partitions_with_coordinates(%partition_coordinates_for);

	return $partition_name;
}   # ----------  end of subroutine get_partitions  ----------


#######################################
# checking if the partition statement
# exists in the nexus file.
#######################################

sub are_partition_sets {
    my ($partition_name) = $file_content =~ m{^\s*set\s+partition\s*=\s*
                                                ( [^;\s]+ )  #Capturing partition name 
                                                \s*;         #The end of a nexus command ';'
                                             }imx;
    return $partition_name;
}   # ----------  end of subroutine are_partitions_set  ----------


###########################################
# getting the partition names and how many
###########################################

sub get_partition_names {
	my ($partition_name) = @_;
	my ($partition_number, $partition_names) = $begin_set =~ m{^\s*partition\s+$partition_name\s*=\s*
                                                                  ( \d+ )     #partition number
                                                                  \s*:\s*
                                                                  ( [^;]+ )   #partition names with commas
                                                                  \s*;        #end of nexus command ';'
                                                               }imx;
    my @partition_names = split /\s*,\s*/mx, $partition_names if ($partition_number and $partition_names);
	
    return $partition_number, @partition_names;
}   # ----------  end of subroutine get_partition_names  ----------


###########################################
# getting partition coordinates for each
# detected and named partition
###########################################

#This is a fast subroutine. For 2000 partitions takes less than 1 sec
sub get_partition_coordinates {
    my (@partition_names) = @_;
    my %partition_coordinates_for;
    my @sites;
    
    foreach my $partition_name (@partition_names) {
        $begin_set =~ /charset\s+$partition_name\s*=\s* ([^;]+) ;/imx;
        
        if ($1) {
            my $charset = $1;
            #if charset contains other than numbers, spaces, dashes, points or \
            if ($charset =~ /([^\s\d\.\\\-])/xm) {
                return 0, "Unexpected charset character [$charset] in [$partition_name] charset";
            }
            else {
                #removing white spaces from set of coordinates
                $charset =~ s/\s*-\s*/-/gxm;
                $charset =~ s/\s+\\/\\/gxm;
                $charset =~ s/\n+/ /gxm; #not necessary
                @sites = get_charset_sites($charset);
            }
        }
        $partition_coordinates_for{$partition_name} = [@sites];
	}
    
    print "\n<!-- getting partitions -->";
    return 1, 'OK', \%partition_coordinates_for;

}   # ----------  end of subroutine get_partition_coordinates  ----------

###########################################
# getting charset sites 
# detected and named partition
###########################################

sub get_charset_sites {
    my ($charset)        = @_;
    my $RECORD_SEPARATOR = q{\s+};
    my $DASH_SEPARATOR   = q{-};
    my @sites;

    my @charsets = split $RECORD_SEPARATOR, $charset;
    foreach my $sites (@charsets) {
        #charset nt1 = 1-.\3;
        if ($sites =~ /(\d+) $DASH_SEPARATOR (\.|\d+) \\ (\d+)/xm) {
            my $ini_interval = $1;
            my $end_interval = $2 eq q{.}? $align_length : $2;
            my $increment    = $3;
            while ($ini_interval <= $end_interval) {
                push @sites, [$ini_interval, $ini_interval];
                $ini_interval += $increment;
            }
        }
        # 628-633 1609-1614 2083-2088 3076-3081 3424-3429 
        elsif ($sites =~ /$DASH_SEPARATOR/xm) {
            while ($sites =~ /(\d+) $DASH_SEPARATOR (\d+) /gxm) {
                my $ini_interval = $1;
                my $end_interval = $2;
                push @sites, [$ini_interval, $end_interval];
            }
        }
        elsif ($sites =~ /^(\d+)$/xm) {
                my $ini_interval = $1;
                my $end_interval = $1;
                push @sites, [$ini_interval, $end_interval];
        }
    }

    #each array element contains initial and final coordinates
    #it can be initial=end
    return @sites;

}   # ----------  end of subroutine get_charset_sites  ----------



##################################################
# creating the partition - i.e. slicing
# the alignment using the coordinates - in
# the desired format
##################################################
# This is the bottleneck. 366" for 2000 partitions
sub create_partitions_with_coordinates {
    my %partition_coordinates_for = @_;
    
  	foreach my $partition_name (keys %partition_coordinates_for) {
        my @sliced_align_objects;
        my @slices = @{$partition_coordinates_for{$partition_name}};

        my $align_out = Bio::AlignIO->new( -file        => ">$path_to_user_folder/$partition_name.phylip",
                                           -format      => "$output_format",
                                           -interleaved => 0,
                                           -line_length => 10,
                                         );
        
        foreach my $slice (@slices) {
            #if true (1) will keep gap-only columns in the newly created slice. All kind gaps.
            my $sliced_align_object = $align_object->slice($slice->[0],
                                                           $slice->[1],
                                                           1);
            push @sliced_align_objects, $sliced_align_object;
        }
        
        if (scalar @sliced_align_objects == 1) {
            $align_out->write_aln($sliced_align_objects[0]);
        }
        else {
            my $concatenated_align_object = concatenate_slices(@sliced_align_objects);
            $align_out->write_aln($concatenated_align_object);
        }
        
        #Workaround to avoid timeout in the server
        print "\n<!-- partition $partition_name created -->";
    }
    
    return ;
}   # ----------  end of subroutine create_partitions_with_coordinates  ----------


###########################################
# concatenating different slices of an 
# alignment
###########################################

sub concatenate_slices {
    my @sliced_align_objects      = @_;
    my $concatenated_align_object = Bio::SimpleAlign->new;

    #getting sequences ids
    my @seq_ids;
    foreach my $seq_object ( $sliced_align_objects[0]->each_seq() ) {
        push @seq_ids, $seq_object->id();
    }

    foreach my $id (@seq_ids) {
        my $new_seq;
        foreach my $align_object (@sliced_align_objects) {
            for ( $align_object->each_seq_with_id($id) ) {
                $new_seq .= $_->seq();
            }
        }
        my $lennewstring = length $new_seq;
        # print "NEW SEQUENCE IS $newstring\n with $lennewstring AA\n";
        my $new_seq_object = new Bio::LocatableSeq(
                                -seq      => "$new_seq",
                                -id       => "$id",
                                -start    => 1,
                                -end      => "$lennewstring",
                                #-alphabet => "$type"
                              );
        
        $concatenated_align_object->add_seq($new_seq_object);
    }
    
    return $concatenated_align_object;
}

1;

