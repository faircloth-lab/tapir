use strict;
use warnings;
use Carp;
use File::Basename qw(basename);

my $tree_depth = 100; #value used as max to match profile and tree
my $correction_factor; #Value used to correct branch lengths 
                       #and now has to be applied to the rates

my $tree_file_name = glob "Tree_*.newick";
if ($tree_file_name) {
    ($correction_factor, $tree_depth) = $tree_file_name =~ /Tree_(\d+)_([\d\.e\+\-]+)\.newick$/xm;
}

my @hyphy_outfiles = glob "*.hyphy_output";
my $allrates_file = "Allrates_file.allrates";
my %rates_summary_for;

open my $ALLRATES, '+>', "$allrates_file"
    or croak " $0 : failed to open input file '$allrates_file' : $!\n";

foreach my $hyphy_outfile (@hyphy_outfiles) {
    my $loci_name = basename $hyphy_outfile;
    $loci_name =~ s/\.hyphy_output$//xm;

    my ($rates_ref, $num_undefined)
            = _access_hyphy_results($hyphy_outfile, $correction_factor);
    my @rates = @{$rates_ref};
    
    my $total_rates = scalar @rates;
    my $total_sites = $total_rates + $num_undefined;

    $rates_summary_for{$loci_name} = [ "$loci_name.hyphy_output",
                                        $total_sites,
                                        $total_rates,
                                        $num_undefined,
                                     ];

    print {$ALLRATES} $loci_name, q{:}, (join q{,}, @rates), "\n";
}

close $ALLRATES
    or carp "$0 : failed to close output file '$allrates_file' : $!\n";

#Creating table with results              
my $table = build_HTML_table(%rates_summary_for);

sub build_HTML_table {
    my ( %rates_summary_for )= @_;
    my $table;
    my $count = 1;

    foreach my $partition (keys %rates_summary_for) {
        $table .=  <<"END_TABLE";
 
                     <tr><td class="col1">$count</td>
                     <td class="col2">$partition</td>
                     <td class="col3">
                     <a href="$rates_summary_for{$partition}[0]">hyphy</a></td>
                     <td class="col4">$rates_summary_for{$partition}[1]</td>
                     <td class="col5">$rates_summary_for{$partition}[2]</td>
                     <td class="col6">$rates_summary_for{$partition}[3]</td>
                     <td class="col7"><input class="part_check" name="_getPI_$partition" type="checkbox" /></td>
                     <td NOWRAP><input class="colorPicker" name="_c_$partition" size="7" value="" type="text" /></td>

END_TABLE

        $count++;
    }

    return $table;
}# ----------  end of subroutine build_HTML_table  ----------


#########################################
# Reading the hyphy file and getting 
# the rates
#######################################

sub _access_hyphy_results {
    my ($file, $correction_factor) = @_;
    my @rates;
    my $num_undefined = 0;

    open my $INFILE, '<', "$file"
        or croak " $0 : failed to open '$file' : $!\n";

    while (my $line = <$INFILE>) {
        next if $line !~ /^\s*Site/xm; 
        #A line with rates from hyphy output Site / Rate
        #Site    1 Total subst =  1.3710 subst, Rate =  0.1523 subst/time, Log(L) -4.0538
        if ($line =~ /Rate\s+=\s+ ( [\d\.]+ )/xm) {
           my $corrected_rate = $1/$correction_factor;
           push @rates, $corrected_rate;
        }
        else {
            ++$num_undefined;
        }
    }
    close $INFILE
        or carp "$0 : failed to close input file '$file' : $!\n";

    return (\@rates, $num_undefined);
}