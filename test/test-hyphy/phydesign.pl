use warnings;
use strict;
use English qw( -no_match_vars );
use Carp;
use Bio::AlignIO;
use Bio::TreeIO;
use IO::String;
use List::Util qw( first );
use Nexus_partitioning;

my $tree_object;
my %taxa_tree;
my $format;
my $root_object;
my $tree_length;
my $one_leaf_node;
my $taxa_seperator = '*=*';
my $sid="bcf";

my $y = format_checking_for_Alignment("../test-data/chr1_918.nex");

eval {
    $tree_object = new Bio::TreeIO(-file => "../test-data/Euteleost.tree",
                                   -format => 'newick',)->next_tree();
    foreach my $leaf ($tree_object->get_leaf_nodes) {
        $taxa_tree{$leaf->id}++; #get leaf name
        $one_leaf_node = $leaf;
    }
    $root_object = $tree_object->get_root_node;
    $tree_length = $tree_object->total_branch_length;
};

my %id_align_for;
%id_align_for = get_taxonID_from_file("IDs_table_$sid.txt");

#my %taxonID_for = create_taxonID_file(sort @taxa);
#my %id_align_for = get_taxonID_from_file("IDs_table_$sid.txt");

#Getting the tree in a string
my $tree_string;
my $io_string = new IO::String(\$tree_string);
my $out = new Bio::TreeIO(-fh     => $io_string,
                          -format => 'newick');
$out->write_tree($tree_object);

#Remove from tree labels or bootstrap values
#i.e., any kind of character between ) and : or ;
$tree_string =~ s/ \) [^:]+ (:|;) /)$1/gxm;

#Printing tree in SVG format
#tree_string_to_svg($tree_string, "Tree.svg");

#getting the depth of the tree
my $tree_depth = $one_leaf_node->depth;

#Rounding number to 3 decimal places or exponencial annotation
$tree_depth = $tree_depth >= 0.001 ? sprintf '%.3f', $tree_depth
                                   : sprintf '%.2e', $tree_depth;

print "Tree depth:\t\t", $tree_depth;
print "\n";
#Branches should be 10(1) or lower
#Total of branches is 2n-3 where n is leaves.
my $total_branches = 2 * (scalar keys %taxa_tree) - 3;
my $mean_branch_length =  $tree_length/$total_branches;
print "Mean branch length:\t", $mean_branch_length;
print "\n";
my $string_length = length int($mean_branch_length + .5);
print "String length:\t\t", $string_length;
print "\n";
my $correction_factor = $string_length > 1 ? ('1' . '0' x $string_length) : 1;
print "Correction factor:\t", $correction_factor;
print "\n";

#dividing each branch by the correction factor
$tree_string =~ s{(?<=:)([\d\.]+)}{$1/$correction_factor}gexm;

#renaming tree leaves with the new tags
#if alignment was uploaded
$tree_string = rename_tree_leaves($tree_string, \%id_align_for);

#writting tree in a file
#with correction factor and tree depth in the name
my $OUTTREE_file_name = "Tree_$correction_factor" . "_$tree_depth.newick";
open my $OUTTREE, '>', "$OUTTREE_file_name"
    or croak  "$0 : failed to open  output file '$OUTTREE_file_name' : $!\n";
print {$OUTTREE} "$tree_string";
close $OUTTREE
    or carp "$0 : failed to close output file '$OUTTREE_file_name' : $!\n";

####################################################
# Renaming the tree leaves using the taxon id from
# the file created in the server.
# To use the treestring ...
# Maybe it is not he best way?
####################################################

sub rename_tree_leaves {
    my ($tree_string, $id_align_for_ref) = @_;
    my %id_align_for = %{$id_align_for_ref};

    foreach my $taxon (keys %id_align_for) {
        $tree_string =~ s/$taxon\b/$id_align_for{$taxon}/xm;
    }

    return $tree_string;
} # ----------  end of subroutine rename_tree_leaves  ----------

####################################################
# Get taxa names from file in the server
####################################################
sub get_taxonID_from_file {
    my ($infile_name) = @_;
    my %id_align_for;

    open  my $INFILE, '<', $infile_name
        or croak  "$0: failed to open  input file '$infile_name' : $!\n";

    while (<$INFILE>) {
        chomp;
        my ($taxon_name, $taxon_ID) = split /\Q$taxa_seperator/xm, $_;
        $id_align_for{$taxon_name} = $taxon_ID;
    }

    close  $INFILE
        or carp "$0: failed to close input file '$infile_name' : $!\n";


    return %id_align_for;
}    # ----------  end of subroutine  ----------

#############################################################################
# Create a file with an ID table. This table will be used
# to rename the sequences for our own names/ids (numbers).
#############################################################################
sub create_taxonID_file {
    my @taxa = @_;
    my %taxonID_for;

    open my $IDTABLE, '>', "IDs_table_$sid.txt"
        or croak "$0 : failed to open  output file 'IDs_table_$sid.txt' : $!\n";

    foreach my $x (0..$#taxa) {
        my $taxon = $taxa[$x];
        my $new_ID =  'X_' . ($x+1) . '_X';
        $taxonID_for{$taxon} = $new_ID;
        print {$IDTABLE} $taxon, $taxa_seperator, $new_ID, "\n";
    }

    close $IDTABLE
        or carp "$0 : failed to close output file 'IDs_table_$sid.txt' : $!\n";

    return %taxonID_for;
}    # ----------  end of subroutine create_taxonID_file  ----------

#####################################
# Checking format from alignment file
#####################################
sub format_checking_for_Alignment {
    my ($filename)= @_;
    my $align_object;
    my @seq_objects;

    #return error if we don't guess any of the 3 formats: Nexus, phylip or fasta.
    my ($align_format, $datatype) = guess_align_format("$filename");
    if ($align_format eq 'unknown') {
        return 0, 'Unknown format. Only NEXUS, PHYLIP and FASTA accepted.';
    }

    my $return_msg = "Alignment with format [$align_format] uploaded.";

    #If format nexus and morphology characters
    if ($align_format eq 'nexus' and $datatype =~ /standard/ixm) {
        my ($ok, $msg) = get_morphology_matrix("$filename");
        return $ok, $msg;
    }

    #Getting alignment object using Bioperl module
    eval{
        #To avoid to get stuck in the creation of the object the alarm is set
        #This can happen with some error in nexus format
        local $SIG{ALRM} = sub { croak "alarm\n" };       # NB \n required
        alarm 60; #Setting alarm to a minute
        $align_object = Bio::AlignIO->new(-file => "$filename",
                                          -format => "$align_format")->next_aln;
        alarm 0; #Removing alarm if object is created without problem
        @seq_objects = $align_object->each_seq();
       };

    #Getting error message from Bioperl or creating one
    #Specially when bioperl doesnt create object neither error
    #Returm errpr amd msg
    if ($EVAL_ERROR) {
        my $error_from_bioperl = $EVAL_ERROR =~ /MSG:\s([^\n]+)\n/xm  ? "$1"
                               : defined $align_object                ? "Not a valid [$align_format] format."
                               : $EVAL_ERROR =~ /alarm/xm             ? "Couldn't process [$align_format] format. Check formating."
                               :                                        'New error'
                               ;

        #It cannot die here so it is carp
        carp "ERROR FROM BIOPERL:: $error_from_bioperl";
        return 0, $error_from_bioperl;
    }

    #Create a file with the an ID table. This table will be used
    #to rename the sequences for our own names/ids (numbers).
    #first getting sequences id
    my $renamed_align_object = Bio::SimpleAlign->new; #used to create the new alignment
    open my $IDTABLE, '>', "IDs_table_$sid.txt"
        or croak "$0 : failed to open  output file 'IDs_table_$sid.txt' : $!\n";

    my %taxa_from_seq; #used to test if seq id is duplicated
    foreach my $x (0..$#seq_objects) {

        my $taxon_name = $seq_objects[$x]->id();
        my $new_seq_id =  'X_' . ($x+1) . '_X';

        print {$IDTABLE} $taxon_name, $taxa_seperator, $new_seq_id, "\n";

        #test for duplicated seq ids
        if (exists $taxa_from_seq{$taxon_name}) {
            return 0, "ID [$taxon_name] is duplicated. It may not be the only one." ;
        }
        $taxa_from_seq{$taxon_name}++;

        #Setting up the new sequence ID
        $seq_objects[$x]->display_id($new_seq_id);
        $renamed_align_object->add_seq($seq_objects[$x]);
    }

    close $IDTABLE
        or carp "$0 : failed to close output file 'IDs_table_$sid.txt' : $!\n";


    #writting out the alignment in phylip format with the renamed sequences
    #if format is nexus, check if there is partition definitions
    if ( $align_format eq 'phylip' || $align_format eq 'fasta' ) {
        create_output_alignment($renamed_align_object, 'phylip');
    }
    elsif ($align_format eq 'nexus') {

        my ($are_partitions, @error_found) = Nexus_partitioning::get_partitions($filename, $renamed_align_object, 'phylip');

        if ($error_found[0]) {
            return 0, $error_found[1];
        }
        elsif (!$are_partitions) {
            #If no partition print whole alignment
            create_output_alignment($renamed_align_object, 'phylip');
            $return_msg .= ' No partitions has been created.'
        }
        else {
            $return_msg .= " Partitions named [$are_partitions] has been created.";
        }


    }

    return 1, $return_msg;
}    # ----------  end of subroutine format_checking_for_Alignment  ----------

####################################################
# Guessing alignment format based on the first line
####################################################
sub guess_align_format {
    my ($infile_name) = @_;

    open  my $INFILE, '<', $infile_name
        or croak  "$0 : failed to open  input file '$infile_name' : $!\n";
    my @file_content = <$INFILE>;
    close  $INFILE
        or carp "$0 : failed to close input file '$infile_name' : $!\n";

    my $first_noempty_line = first {$_ =~ /\S/xm} @file_content;
    #We only get defined values for the next two variables if file equal
    #nexus format
    my $datatype_line = first {$_ =~ /datatype/ixm} @file_content;
    my $datatype = 'unidentified';
    if ($datatype_line) {
        ($datatype) = $datatype_line =~ /datatype\s*=\s*(\S+)/ixm;
    }

    my $align_format = $first_noempty_line =~ /^\s*>/xm              ? 'fasta'
                     : $first_noempty_line =~ /^\s*\d+\s+\d+\s*$/xm  ? 'phylip'
                     : $first_noempty_line =~ /^\s*\#nexus/ixm       ? 'nexus'
                     :                                                 'unknown'
                     ;

    return $align_format, $datatype;
}   # ----------  end of subroutine guess_align_format  ----------

####################################################
# Writting in user folder the new alignment
####################################################
sub create_output_alignment {
    my ($align_object, $format) = @_;
    my $align_out = Bio::AlignIO->new(-file   => ">Alignment.$format" ,
                                      -format => "$format",
                                      -interleaved => 0,
                                      -line_length  => 10,);
    $align_out->write_aln($align_object);

    return ;
}   # ----------  end of subroutine create_output_alignment  ----------