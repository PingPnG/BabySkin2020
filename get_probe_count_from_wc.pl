#! /usr/bin/perl -w
##########################################################################
#### Author: Ping HU
#### 
##########################################################################
use strict;
my $usage="$0  [file] first line is title\n"; 
die $usage unless @ARGV ==1;
my ($file) = @ARGV;

open (A, "<$file")||die "could not open $file\n";

while (my $a=<A>){
    chomp $a;
    my @tmp=split /\s+/, $a;
    print $tmp[4], "\t", $tmp[1]-1, "\n";
}
close A;
