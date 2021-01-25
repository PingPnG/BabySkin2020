#! /usr/bin/perl -w
##########################################################################
#### Author: Ping HU
#### 
##########################################################################
use strict;
my $usage="$0  [file][kofile] first line is title\n"; 
die $usage unless @ARGV ==2;
my ($file, $kofile) = @ARGV;

open (A, "<$kofile")||die "could not open $file\n";
my %ko;
while (my $a=<A>){
    chomp $a;
    my @tmp=split /\t/, $a;
    if($tmp[1] =~ /ko:(\S+)/){
	$tmp[1]=$1;
    }
    $ko{$tmp[0]}=$tmp[1];
}
close A;


open (A, "<$file")||die "could not open $file\n";
my $a=<A>;
chomp $a;
for($a){s/\r//gi;}
my $tt=split /\t/, $a;
print "KOID\tPvalue\tFC\tgene\n";
while (my $a=<A>){
    chomp $a;
    my @tmp=split /\t/, $a;
    if(defined $ko{$tmp[0]}){
	print $ko{$tmp[0]},"\t",  $tmp[1], "\t", $tmp[2], "\t", $tmp[0], "\n";
    }else{
	print $tmp[0],"\t",  $tmp[1], "\t", $tmp[2], "\t", $tmp[0], "\n";
    }
}
close A;

