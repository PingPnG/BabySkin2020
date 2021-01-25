#! /usr/bin/perl -w
##########################################################################
#### Author: Ping HU
#### 
##########################################################################
use strict;
my $usage="$0  [file] first line is title\n"; 
die $usage unless @ARGV ==1;
my ($file) = @ARGV;
my %data;
open (A, "<$file")||die "could not open $file\n";

while (my $a=<A>){
    chomp $a;
    my @tmp=split /\#/, $a;
    foreach my $i (@tmp){
	$data{$i}=1;
    }
}
close A;

foreach my $i (keys %data){
    print $i, "\n";
}
