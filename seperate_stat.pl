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
my $x=<A>;
for($x){s/\r//gi;}
chomp $x;
my @tmp=split /\t/, $x;
close A;
my %match;
for (my $i=1; $i<=@tmp; $i++){
    if($tmp[$i] =~ /pt\_(\S+)/){
	$match{$1}{"p"}=$i+1;
    }elsif($tmp[$i] =~ /TFC\_(\S+)/){
	$match{$1}{"fc"}=$i+1;
    }
}
my $nn="cut -f1,";
foreach my $i(sort keys %match){
    my $n=0;
    if(defined $match{$i}{"p"}){
	print "cut -f1,", $match{$i}{"p"}, " $file>tmp.$file.$i.p\n";
	$n++;
    }
    if(defined $match{$i}{"fc"}){
	print "cut -f", $match{$i}{"fc"}, " $file>tmp.$file.$i.fc\n";
	$n++;
    }
    if($n==2){
	print "paste tmp.$file.$i.p tmp.$file.$i.fc >$file.$i.stat\n";
    }
}

print "rm -rf tmp.$file.*\n";
print "mkdir ./stat\n";
print "mv $file.*.stat ./stat/\n";
