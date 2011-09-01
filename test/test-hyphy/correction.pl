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