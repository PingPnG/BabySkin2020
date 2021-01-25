#! /usr/bin/perl -w
##########################################################################
#### Author: Ping HU
#### 
##########################################################################
use strict;
my $long="/home/ping/db/affy/NA36/HG-U219.na36.annot.csv.long";
my $syn="/home/ping/db/affy/NA36/HG-U219.na36.annot.csv.syn";
my $loc="/home/ping/db/affy/NA36/HG-U219.na36.annot.csv.loc";
####forallgene
#my $preEpi="/home/ping/project/Oral2527_2537/GSS2527_2537/src/GSS2527.pre30";
#my $preDerm="/home/ping/project/Oral2527_2537/GSS2527_2537/src/GSS2537.pre30";
####forsinglegene
#my $preEpi="/home/ping/project/Oral2527_2537/GSS2527_2537/signal_file/Epi.1gene.30.probe";
#my $preDerm="/home/ping/project/Oral2527_2537/GSS2527_2537/signal_file/Derm.1gene.30.probe";
my $preEpi="../src/allprobe";
my $ko="/home/ping/db/affy/NA36/HG-U219.na36.ko";
my $go="/home/ping/db/affy/NA36/HG_U219.na36.go";
my ($file, $name) = @ARGV;
open (B, "<$file")|| die "could not open $file\n";
#open (C, ">$name.heatjob")|| die "could not open $name.heatjob\n";
#open (D, ">$name.enrichjob")||die "could not open $name.enrichjob\n";
#open (F, ">$name.pathwayjob")|| die "could not open $name.pathwayjob\n";
open (E, ">$name.keggjob")|| die "could not open $name.keggjob\n";
while (my $line =<B>){
    chomp $line;
    my $pre=$preEpi;
    #   if($line =~ /derm/gi){
    #	$pre=$preDerm;
    #}
    # print F "perl /home/ping/project/Oral2527_2537/GSS2527_2537/src/locate_pathway3_filter_new.pl  ../stat/$line $syn  $long $pre\n";

# print D "perl /home/ping/project/Oral2527_2537/GSS2527_2537/src/match_KOID.pl  ../stat/$line $ko >../stat/$line.ko\n";
# print D "Rscript /home/ping/project/Oral2527_2537/GSS2527_2537/src/enrichment_test.R ../stat/$line.ko\n";
#print E "perl /home/ping/project/Oral2527_2537/GSS2527_2537/src/test_kegg_compound_new.pl $ko ../stat/$line  $syn $long  $pre \n";
    #   print C "perl /home/ping/project/Oral2527_2537/GSS2527_2537/src/write_html_new_cor.pl $pre $go 0.1 ../stat/$line $long $syn -l $loc $pre \n";
 
  print E "perl ../src/color_map_box_ko_filter_cor.pl $ko ../stat/$line $syn  $long -l $loc $pre \n";


}




close B;
#close C;

#close D;
#close F;
close E;

system ("chmod 777 *job\n");
