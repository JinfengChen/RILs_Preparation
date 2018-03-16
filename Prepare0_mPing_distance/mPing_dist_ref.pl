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
my %type;
open IN, "$file" or die "$!";
while(<IN>){
    chomp $_;
    next if ($_=~/^$/);
    my @unit=split("\t",$_);
    push @{$hash{$unit[0]}}, $unit[3];
    my $flag = 'Non_Ref';
    if ($unit[8]=~/Reference\-only/){
       $flag = "ref_only";
    }elsif($unit[8]=~/Shared/){
       $flag = "shared";
    }
    $type{"$unit[0]\_$unit[3]"} = $flag;
}
close IN;

foreach my $c (sort keys %hash){
    my @pos = sort {$a <=> $b} @{$hash{$c}};
    my $dist_first = $pos[1] - $pos[0];
    my $index = "$c\_$pos[0]";
    print "$c:$pos[0]\t$dist_first\t$type{$index}\n";
    for (my $i=1; $i<@pos-1; $i++){
        my $dist1 = $pos[$i+1] - $pos[$i];
        my $dist2 = $pos[$i] - $pos[$i-1];
        my $dist  = $dist1 < $dist2 ? $dist1 : $dist2;
        my $index = "$c\_$pos[$i]";
        print "$c:$pos[$i]\t$dist\t$type{$index}\n"; 
    }
    my $dist_last = $pos[$#pos] - $pos[$#pos-1];
    my $index = "$c\_$pos[$#pos]";
    print "$c:$pos[$#pos]\t$dist_last\t$type{$index}\n"; 
}
}
