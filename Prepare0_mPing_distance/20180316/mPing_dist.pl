#!/usr/bin/perl
use Getopt::Long;
use File::Basename;
use Data::Dumper;
use FindBin qw($Bin);


GetOptions (\%opt,"input:s","help");


my $help=<<USAGE;

perl mPing_dist.pl --input HEG4_2.3.mPing.20X.mping.all_inserts.gff > mPing_dist.txt

USAGE


if ($opt{help} or keys %opt < 1){
    print "$help\n";
    exit();
}

readtable($opt{input});

sub readtable
{
my ($file)=@_;
my %hash;
my %strand;
open IN, "$file" or die "$!";
while(<IN>){
    chomp $_;
    next if ($_=~/^$/);
    my @unit=split("\t",$_);
    push @{$hash{$unit[0]}}, $unit[3];
    $strand{"$unit[0]\.$unit[3]"} = $unit[5];
}
close IN;

my $cutoff = 50000000;
#my $cutoff = 100000;
open OUT1, ">mPing_dist2.txt" or die "$!";
open OUT2, ">mPing_dist2.50Mb.list" or die "$!";
#open OUT2, ">mPing_dist2.100kb.list" or die "$!";
#open OUT1, ">mPing_dist_RIL_AF0.1.txt" or die "$!";
#open OUT2, ">mPing_dist_RIL_AF0.1.50Mb.list" or die "$!";
foreach my $c (keys %hash){
    my @pos = sort {$a <=> $b} @{$hash{$c}};
    my $dist_first = $pos[1] - $pos[0];
    print OUT1 "$c\t$pos[0]\t$dist_first\n";
    if ($dist_first <= $cutoff){
        my $mping1 = "$c\.$pos[0]";
        my $mping2 = "$c\.$pos[1]";
        print OUT2 "$mping1\t$mping2\t$dist_first\t$strand{$mping1}\t$strand{$mping2}\n";
    }
    for (my $i=1; $i<@pos-1; $i++){
        my $dist1 = $pos[$i+1] - $pos[$i];
        my $dist2 = $pos[$i] - $pos[$i-1];
        if ($dist1 < $dist2){
            my $dist = $dist1;
            print OUT1 "$c\t$pos[$i]\t$dist\n";
            if ($dist <= $cutoff){
                my $mping1 = "$c\.$pos[$i]";
                my $mping2 = "$c\.$pos[$i+1]";
                print OUT2 "$mping1\t$mping2\t$dist\t$strand{$mping1}\t$strand{$mping2}\n";
            }
        #}elsif($dist1 > $dist2){
        }else{
            my $dist = $dist2;
            print OUT1 "$c\t$pos[$i]\t$dist\n";
            if ($dist <= $cutoff){
                my $mping1 = "$c\.$pos[$i-1]";
                my $mping2 = "$c\.$pos[$i]";
                print OUT2 "$mping1\t$mping2\t$dist\t$strand{$mping1}\t$strand{$mping2}\n";
            }
        }
        #my $dist  = $dist1 < $dist2 ? $dist1 : $dist2;
        #print OUT1 "$c\t$pos[$i]\t$dist\n";
        #if ($dist <= 50000000){
        #    my $mping1 = "$c\.$pos[$i]";
        #    my $mping2 = "$c\.$pos[$i+1]";
        #    print OUT2 "$mping1\t$mping2\t$dist\t$strand{$mping1}\t$strand{$mping2}\n";
        #} 
    }
    my $dist_last = $pos[@pos-1] - $pos[@pos-2];
    print OUT1 "$c\t$pos[@pos-1]\t$dist_last\n";
    if ($dist_last <= $cutoff){
        my $mping1 = "$c\.$pos[@pos-2]";
        my $mping2 = "$c\.$pos[@pos-1]";
        print OUT2 "$mping1\t$mping2\t$dist_last\t$strand{$mping1}\t$strand{$mping2}\n";
    }
}
close OUT1;
close OUT2;
}
