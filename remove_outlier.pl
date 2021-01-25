#! /usr/bin/perl -w
##########################################################################
#### Author: Ping HU
#### 
##########################################################################
use strict;
my $usage="$0  [outlier_id_file][signal_file] first line is title\n"; 
die $usage unless @ARGV ==2;
my ($outfile, $file) = @ARGV;
#my $outfile="/home/ping/project/HairWinterClinical/XingtaoAnalysis/BC4.outlier.txt";
my %id;
open (A, "<$outfile")||die "could not open $outfile\n";
my $x=<A>;
while (my $a=<A>){
    chomp $a;
    my @tmp=split /\t/, $a;
    $id{$tmp[0]}=1;
}
close A;

#my $file="GSS1381.matchdata";
open (A, "<$file")||die "could not open $file\n";
my $x=<A>;
if($x  =~ "Constructed from biom file"){
    $x=<A>;
}
chomp $x;
for($x){s/\r//gi;}
my @tmp=split /\t/, $x;
my $command= "cut -f1";
for(my $i=1; $i<@tmp; $i++){
   
    if(defined $id{$tmp[$i]}){
	next;
    }else{
	my $y=$i+1;
	$command.= ",". $y;
    }
}
system("$command $file > $file.no_outlier\n");
close A;
