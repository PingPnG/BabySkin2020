#!/usr/bin/perl -w
use strict;
my $usage="$0  [file] first line is title\n"; 
die $usage unless @ARGV ==1;
my ($file) = @ARGV;
open (Z, ">$file.truefc")||die "could not open $file.truefc\n";
open (X, "<$file")||die "could not open $file\n";
my $title=<X>;
chomp $title;
for($title){s/\r//gi;}
my @tmp=split /\t/, $title;
print Z "Gene\t", $tmp[1], "\tTrueFC.from.", $tmp[2], "\n";
while (my $line = <X>){
    chomp $line;
    for($line){s/\r//gi;}
    my @tmp=split /\t/, $line;
    print Z $tmp[0], "\t",$tmp[1], "\t";
    my $x="NA";
    if ($tmp[2 ]>=0){
	$x=2**$tmp[2];
    }else{
	$x=-2**(-$tmp[2 ]);
    }
    print Z $x, "\n";
   
}
close X;
close Z;
