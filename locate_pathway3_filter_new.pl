#! /usr/bin/perl -w
############Ping Hu#####################
#
use strict;
use Getopt::Long;


my $PLOT_LINK2="http://mvic-biotech.na.pg.com/projects/Hair_Growth/MDEScalpHFGraph/";

my $PLOT_LINK=$PLOT_LINK2;


my ($opt_l, $opt_g);
GetOptions(
	   "l=s" => \$opt_l,
	   "g" => \$opt_g
	   );

my $usage = "$0 stat_file Acronyme_file long_file

 (1) output file: Experiment_name.html and Experiment_name
 (2) stat file format: affy\tP_value\tFold_change
 (3) Acronyme file format: affy\tAcronym
 (4) Long annotation file name: affy\tlong
  
  -l chromosome mapping location file

\n";

die $usage unless @ARGV == 4;
my ($statfile, $acr_file, $long_file, $filter) = @ARGV;
my $log;
if ($opt_g){
    $log=1;
}

my $outname= $statfile;
if ($statfile =~  /([^\/]+)\s*$/ ){
    $outname=$1;
}
if ($outname =~/(\S+)\.stat/){
    $outname=$1;
}
#####################
# constant
####################
my $COLWID=180;
#my $pathfile="/home/ping/db/pathway/data/pathway-mammal.list";
#my $pathfile="/home/ping/db/pathway/data/pathway-mammal-nonKEGG.list";
my $pathfile="/home/ping/code/all_pathway.map";
#my $long_file="/home/ping/db/affy_chip/Dog.long.txt";
#my $acr_file="/home/ping/db/affy_chip/Dog.FA.txt";
#######################
# step 1. Load data
#######################
my %location;
if ($opt_l){
    open (MMMM, "<$opt_l")or die "could not open $opt_l\n ";
    while (my $line=<MMMM>){
	chomp $line;
	my @tmp=split /\t/, $line;
	$location{$tmp[0]}=$tmp[1];
	   # print STDERR $tmp[0], "\t", $tmp[1], "\n";
    }
    close MMMM;
}

my %acr; #from affy to acronym
fill_hash(\%acr, $acr_file);
my %long;
fill_hash(\%long, $long_file);
my %gene; #from acronym to affy id
foreach my $aff (keys %acr){
    my $aa=$acr{$aff};
    #print STDERR $aa, "*****", $acr{$aff}, "\n";
    if (! defined $gene{$aa}){ 
	$gene{$aa}=$aff;            
    }else{
	$gene{$aa}.="\t".$aff; 
    }
}

my (%FFF, %PPP);

if (defined $log){
   fill_stat_log_filter(\%PPP, \%FFF, $statfile, $filter); 
#fill_stat_log(\%FFF, \%PPP, $statfile);
}else{
    fill_stat_filter(\%PPP, \%FFF, $statfile, $filter);
  #fill_stat(\%FFF, \%PPP, $statfile);  
}
my %list_gene; #genes in the data set acr->aff
foreach my $aff(keys %PPP){
    if (defined $acr{$aff}){
	$list_gene{$acr{$aff}}=$aff;
    }
}

my %path;
my %path_name;
my %path_n;
my %chip_n;
my %data_n;
my %path_gene;
my %path_link;
my %path_genelink;
open(IN,"<$pathfile") or die "Couldn't open  $pathfile";
while(my $in = <IN>){
    chomp $in; 
    my @tmp=split /\t/, $in;
    $tmp[3]=lc($tmp[3]);
   # print $tmp[3], "\n";
    if (scalar @tmp <6){next;}
    $path_name{$tmp[0]}=$tmp[2];
    for(my $i=4; $i<@tmp; $i++){
	if ($tmp[$i] eq ""){next;}
	$path_gene{$tmp[0]}{$tmp[$i]}=$tmp[$i];
	$path_genelink{$tmp[0]}{$tmp[$i]}="http://en.wikipedia.org/wiki/".$tmp[$i];
	if (! defined $path_n{$tmp[0]}){
	    $path_n{$tmp[0]}=1;
	}else{
	    $path_n{$tmp[0]}++;
	}
	if (defined $gene{$tmp[$i]}){
	    #print $tmp[3], "\t", $gene{$tmp[3]}, "\n";
	    if (! defined $chip_n{$tmp[0]}){
		$chip_n{$tmp[0]}=1;
	    }else{
		$chip_n{$tmp[0]}++;
	    }
	}
	if (defined $list_gene{$tmp[$i]}){
	    $path{$tmp[0]}{$tmp[$i]}=$list_gene{$tmp[$i]};
	    if (! defined $data_n{$tmp[0]}){
		$data_n{$tmp[0]}=1;
	    }else{
		$data_n{$tmp[0]}++;
	    }
	}
    }
    
    $path_link{$tmp[0]}=$tmp[1];
    
 
  
  
}
close IN;

#######################
# step 2 print html files
######################
my $outfile="$outname.html";
system ("mkdir $outname");
my $header=get_header($outname);
open (B, ">$outfile")or die "could not open $outfile\n ";
print B $header;

my %part; #sort 
foreach my $pp(keys %path){
    my $pname=$path_name{$pp};
    my $header2= $pp."\t".$path_name{$pp}."\n";
    my $N2=$chip_n{$pp};
    my $N3=$data_n{$pp};
    if ($N3<5){next;}
    my $N4=0;
    my $up=0;
    my $down=0;
    my $sigup=0;
    my $sigdown=0;
    foreach my $gg (keys %{$path{$pp}}){
	my @ll=split "\t", $path{$pp}{$gg};
	my $c=0;
	my $gup=0;
	my $gdown=0;
	my $siggup=0;
	my $siggdown=0;
	foreach my $probe (@ll){
	    print STDERR $probe,"#";
	    if (defined $FFF{$probe}){
		for ($PPP{$probe}){
		    s/\</0/gi;
		}
		if ($PPP{$probe}<=0.05){
		    $c=1;
		    if ($FFF{$probe}>1){$siggup++;}
		    elsif($FFF{$probe}<1){$siggdown++;}
		}
		if ($FFF{$probe}>1){$gup++;}
		elsif($FFF{$probe}<1){$gdown++;}
	    }elsif (defined $PPP{$probe}){
		if ($PPP{$probe}<=0.05){$c=1;}
	    }
	}
	$N4+=$c;
	if ($gup>0){$up++;}
	if ($gdown>0){$down++;}	
	if ($siggup>0){$sigup++;}
	if ($siggdown>0){$sigdown++;}
    }
    my $P2= int (100*$N4/$N3);
    my $PU= int (100*$up/$N3);
    my $PD= int (100*$down/$N3);
    my $sigPU=0;
    my $sigPD=0;
    #####add in forced cutoff######
    #if (($N4>=3)&&($P2>=5)){
    if($N4>=1){
	$sigPU= int (100*$sigup/$N4);
	$sigPD= int (100*$sigdown/$N4);
    
	my $sub_file_name="$outname/ITEM_$pp.html";
	write_gene_html($pp,$pname, $outname, $sub_file_name, \%{$path{$pp}}, \%long, $path_link{$pp},\%{$path_genelink{$pp}} );
    
	####here is the score and color code##############
	my $span_class="style=\"background-color:C0C0C0\"";
	my $span_class2="style=\"background-color:C0C0C0\"";
	my $score=0;
	if ($PU >$PD){
	    $span_class2="style=\"background-color:FF9999\"";
	    if ($PU>=90){
		$span_class2="style=\"background-color:FF3333\"";
		$score=$score+2;
	    }elsif($PU>=70){
		$span_class2="style=\"background-color:FF6666\"";
		$score=$score+1;
	    }
	}elsif($PD>$PU){
	    $span_class2="style=\"background-color:CCFFFF\"";
	    if ($PD>=90){
		$span_class2="style=\"background-color:0000FF\"";
		$score=$score+2;
	    }elsif($PD>=70){
		$span_class2="style=\"background-color:0099FF\"";
		$score=$score+1;
	    }	
	}
	
	my $span_class4="style=\"background-color:C0C0C0\"";
	if ($sigPU >$sigPD){
	    $span_class4="style=\"background-color:FF9999\"";
	    if ($sigPU>=90){
		$span_class4="style=\"background-color:FF3333\"";
		$score=$score+2;
	    }elsif($sigPU>=70){
		$span_class4="style=\"background-color:FF6666\"";
		$score=$score+1;
	    }
	}elsif($sigPD>$sigPU){
	    $span_class4="style=\"background-color:CCFFFF\"";
	    if ($sigPD>=90){
		$span_class4="style=\"background-color:0000FF\"";
		$score=$score+2;
	    }elsif($sigPD>=70){
		$span_class4="style=\"background-color:0099FF\"";
		$score=$score+1;
	    }	
	}
	
	$score=$score+$N4;
	if ($P2>=50){
	    if($P2<70){
		$span_class="style=\"background-color:FF3366\"";
		$score=$score+2;
	    }else{
		$span_class="style=\"background-color:FF0000\"";
		$score=$score+3;
	    }
	}
	my $line_content=
	    "<td width=\"$COLWID\"><a href=\"$sub_file_name\">\n$pname</a></td>\n"
	    ."<td><a href=\"".$path_link{$pp}."\">\n$pp</a></td>\n"
	    ."<td $span_class align=\"CENTER\">$P2</td>\n"
	    
	    ."<td $span_class2 align=\"CENTER\">$PU</td>\n"
	    ."<td $span_class2 align=\"CENTER\">$PD</td>\n"
	    ."<td $span_class4 align=\"CENTER\">$sigPU</td>\n"
	    ."<td $span_class4 align=\"CENTER\">$sigPD</td>\n"   
	    ."<td align=\"CENTER\">$N4</td>\n"
	    ."<td align=\"CENTER\">$N2:$N3:$N4:$up:$down:$sigup:$sigdown</td>\n"
	    ."</tr>\n"
	    ."<tr>\n";
	print STDERR  "\t",$outname, "\t",$pp, "\t",  $pname,"\t",$N3,"\t",  $N4, "\t", $sigup, "\t", $sigdown, "\t",  $sub_file_name,   "\n";
	if (! defined $part{$score}){
	    $part{$score}=$line_content;
	}else{
	    $part{$score}.=$line_content;
	}
    }
    
}
foreach my $score (sort {$b<=>$a} keys %part){
    print B $part{$score};
}
print B time_string(), "\n";
print B "</body></HTML>\n";
close B;


sub write_gene_html{
    my ($pid, $pathway_name,$TITLE, $file_name, $path_pp_hash, $long_hash, $link, $path_genelink_hh) =@_;
    open (DDDD, ">$file_name")or die "could not open $file_name\n ";
    my $header=get_small_header($pathway_name, $TITLE, $link);   
    print DDDD $header;
    my $pic=$link;
    if ($pid=~ /KEGG(\d+)/){
	$pic="http://www.genome.ad.jp/kegg/pathway/map/map$1.gif";
	$link="http://www.genome.ad.jp/kegg/pathway/map/map$1.html";
    }elsif ($link=~/biocarta/ ){
	for ($pic){s/\.asp/\.gif/i;}
    }
    print DDDD "<P><IMG src=\"$pic\">\n";

    print DDDD "\n\n<tr><td><font color=\"#000099\">Gene</font></td><td><font color=\"#000099\">AffyProbe</font></td><td><font color=\"#330033\">P</font></td><td><font color=\"#996600\">Fold</font></td><td><font color=\"#003300\">Description</font></td>";
    if ($opt_l){
	print DDDD "<td><font color=\"#006600\">Location</font></td>\n";	
    }
    print DDDD "</tr>\n";
    my %order;
    my $other="";
    foreach my $gene_name (keys %{$path_pp_hash}){
	my @ll=split "\t", $$path_pp_hash{$gene_name};
	my $lname="";
	my $llink="";
       

	if (defined $$path_genelink_hh{$gene_name}){ $llink=$$path_genelink_hh{$gene_name};}	
	foreach my $name(@ll){
	   
	    if (defined $$long_hash{$name}){ $lname=$$long_hash{$name};}
	    my $P="N/A";
	    my $FC="N/A";
	    if (defined $PPP{$name}){$P=$PPP{$name};}
	    if (defined $FFF{$name}){$FC=$FFF{$name};}
	    my $span_class="style=\"background-color:f0f0f0\"";
	    my $span_p="style=\"background-color:C0C0C0\"";
	    if ($FC ne "N/A"){
		if ($FC >4){
		    $span_class="style=\"background-color:f06464\""; 
		}elsif($FC>2){
		    $span_class="style=\"background-color:f09292\"";
		}elsif($FC>1.2){
		    $span_class="style=\"background-color:f0c1c1\"";
		}elsif($FC>0.8){
		    $span_class="style=\"background-color:f0f0f0\"";
		}elsif($FC>0.5){
		    $span_class="style=\"background-color:c1c1f0\"";
		}elsif($FC>0.25){
		    $span_class="style=\"background-color:9292f0\"";
		}else{
		    $span_class="style=\"background-color:6464f0\"";
		}
	    }
	    if ($P ne "N/A"){
		if ($P <0.0001){ 
		    $span_p="style=\"background-color:ff0000\""; 
		}elsif($P<0.001){
		    $span_p="style=\"background-color:ff9700\"";
		}elsif($P<0.01){
		    $span_p="style=\"background-color:ffca00\"";
		}elsif($P<0.1){
		    $span_p="style=\"background-color:fffd00\"";
		}else{
		    $span_p= "style=\"background-color:c0c0c0\"";
		} 
	    }
	    
	    my $gene_line="<tr>\n". 
	    "<td ><font color=\"#006600\"><a href=\"javascript:void(0)\" onclick=\"gsr(\'$gene_name\');\">$gene_name</a></font></td>\n". 
	    "<td ><a href=\"javascript:void(0)\" onclick=\"gar(\'$name\');\">$name</a> </td>\n"."<td $span_p><font color=\"#330033\">$P</font></td>\n".
	    "<td $span_class><font color=\"#996600\">$FC</font></td>\n".
	    "<td><font color=\"#003300\"><a href=\"$llink\">$lname</a></font></td>\n".
	    "</tr>\n";
	    if ($opt_l){
		my $loc="N/A";
		if (defined $location{$name}){
		    $loc= $location{$name}
		} 
		$gene_line="<tr>\n". 
		    "<td ><font color=\"#006600\"><a href=\"javascript:void(0)\" onclick=\"gsr(\'$gene_name\');\">$gene_name</a></font></td>\n". 
		    "<td ><a href=\"javascript:void(0)\" onclick=\"gar(\'$name\');\">$name</a> </td>\n".	                           
		    "<td $span_p><font color=\"#330033\">$P</font></td>\n".
		    "<td $span_class><font color=\"#996600\">$FC</font></td>\n".
		    "<td><font color=\"#003300\"><a href=\"$llink\">$lname</a></font></td>\n<td><font color=\"#006600\">$loc</font></td>\n</tr>\n";
	}

	    if ($FC ne "N/A"){
		if (! defined $order{$FC}){$order{$FC}=$gene_line;}
		else{
		    $order{$FC}.=$gene_line;
		}
	    }else{
		$other.=$gene_line;
	    }

	
	}
    }
    foreach my $i (sort {$a<=>$b} keys  %order){
        print DDDD $order{$i};
    }
    print DDDD $other;
    print DDDD "</table>\n";
    print DDDD time_string(), "\n";
    print DDDD "</BODY></HTML>\n";
    close DDDD;
    return;
}



 
sub fill_hash{
    my ($h, $file)=@_;
    open (X, "$file")||die "could not open $file\n";
    while (my $line = <X>){
	chomp $line;
	my @tmp=split /\t/, $line;
	if (@tmp<2){next;}
	$$h{lc($tmp[0])}=lc($tmp[1]);
    }
    close X;
    return;
}



sub get_header{
    my ($TITLE)=@_;
    my $MAIN_HEADER="<html><head><title>$TITLE</title>
	<script language=\"JavaScript\">
	function gsr( gene ) {
	    link = \'http://en.wikipedia.org/wiki/' + gene;
	    window.open( link, \'popup\', \'width=900,height=700,resizable=yes,scrollbars=yes\' );
	}
	
       function gar( gene ) {
	    link = \'$PLOT_LINK\' + gene + '.png';
	    window.open( link, \'popup\', \'width=900,height=700,resizable=yes,scrollbars=yes\' );
	}
        function fue( gene ) {
	    link = \'$PLOT_LINK2\' + gene + '.png';
	    window.open( link, \'popup\', \'width=900,height=700,resizable=yes,scrollbars=yes\' );
	}

	var newWin;
	function popupWindow( term, geneList, genes ) {
	    if( newWin != null ) {
			newWin.close();
		}
        newWin = window.open(\'http://yahoo.com\',\'popup\',\'width=900,height=700,resizable=yes\'); 
    }
	</script></head><body>

   <h1><center>$TITLE</center></h1>

<table>
<tr><td valign=top><table border=\"1\" >
<tr><td align=\"CENTER\"><b>Column</b></td><td align=\"CENTER\"><b>Column Name</b></td></tr>
<tr><td align=\"LEFT\">N2</td><td align=\"LEFT\">Pathway Gene Number in Annotation Set </td></tr>
<tr><td align=\"LEFT\">N3</td><td align=\"LEFT\">Pathway Gene Number in Data</td></tr>
<tr><td align=\"LEFT\">N4</td><td align=\"LEFT\">Pathway Gene Number with p <= 0.05 </td></tr>
<tr><td align=\"LEFT\">P2%</td><td align=\"LEFT\">Percentage in data N4/N3 %</td></tr>
<tr><td align=\"LEFT\">UP</td><td align=\"LEFT\">Total Pathway Gene Number Upregulated </td></tr>
<tr><td align=\"LEFT\">DN</td><td align=\"LEFT\">Total Pathway Gene Number Downregulated </td></tr><tr><td align=\"LEFT\">UP</td><td align=\"LEFT\">Total Pathway Gene Number Upregulated </td></tr>
<tr><td align=\"LEFT\">DN</td><td align=\"LEFT\">Total Pathway Gene Number Downregulated </td></tr>
/td></tr><tr><td align=\"LEFT\">sigUP</td><td align=\"LEFT\">Total Pathway Gene Number significantly Upregulated p<=0.05 </td></tr>
<tr><td align=\"LEFT\">sigDN</td><td align=\"LEFT\">Total Pathway Gene Number significantly Downregulated p<=0.05 </td></tr>
<tr><td align=\"LEFT\">U%</td><td align=\"LEFT\">Percentage upregulated in data UP/N3 %</td></tr>
<tr><td align=\"LEFT\">D%</td><td align=\"LEFT\">Percentage downregulated in data DOWN/N3 %</td></tr>
<tr><td align=\"LEFT\">sigU%</td><td align=\"LEFT\">Percentage significantly upregulated in data sigUP/N4 %</td></tr>
<tr><td align=\"LEFT\">sigD%</td><td align=\"LEFT\">Percentage significantly downregulated in data sigDOWN/N4 %</td></tr>
</table>
<table border=\"0\" cellspacing=\"0\" cellpadding=\"0\" width=\"100%\">
<tr><td class=\"tdmouseout\">&nbsp;</td>\n</tr></table></td>
<td valign=top><table border=\"1\">
<tr><th>Heat Map Color Legend</th></tr>
<tr><td style=\"background-color:CCCCCC\">P2<50%</td></tr>
<tr><td style=\"background-color:FFfd00\">P2>=50% N4<5</td></tr>
<tr><td style=\"background-color:FF9700\">P2<50% N4>=5</td></tr>
<tr><td style=\"background-color:FF3366\">70%>P2>=50% N4>=5</td></tr>
<tr><td style=\"background-color:FF0000\">P2>=70% N4>=5</td></tr>
</table></td>
<td valign=top>\n<table border=\"1\">
<tr><th>UP-DOWN Color Legend</th></tr>
<tr><td style=\"background-color:FF3333\">UP% >=90%</td></tr>
<tr><td style=\"background-color:FF6666\">70%<= UP%<90%</td></tr>
<tr><td style=\"background-color:FF9999\">UP regulated </td></tr>
<tr><td style=\"background-color:0000FF\">DOWN% >=90% </td></tr>
<tr><td style=\"background-color:0099FF\">90%> DOWN% >=70% </td></tr>
<tr><td style=\"background-color:CCFFFF\">DOWN Regulated</td></tr>
</table></td>
</tr></table><p><p>
<script type=\"text/javascript\" src=\"http://artemis.na.pg.com/sorttable.js\"></script>
<table border=\"1\" align=\"CENTER\" class=\"sortable\" id=\"802043850\">

<tr><td width=$COLWID><b>Pathway Name</b></td><td><b>ID</b></td>
<td><b>P2%</b></td><td><b>U%</b></td><td><b>D%</b></td><td><b>sigU%</b></td><td><b>sigD%</b></td>
<td><b>N4</b></td><td><b>N2:N3:N4:UP:DN:sigUP:sigDOWN</b></td></tr>

<tr>
";
 return $MAIN_HEADER;
}

sub get_small_header{
    my ($pathway_name,$TITLE, $link)=@_;
    my $header="<html><head><title>$pathway_name</title>
	<script language=\"JavaScript\">
	function gsr( gene ) {
	    link = \'http://en.wikipedia.org/wiki/'+ gene;
		window.open( link, \'popup\', \'width=900,height=700,resizable=yes,scrollbars=yes\' );
	}

	function gar( gene ) {
	    link = \'$PLOT_LINK\' + gene + '.png';
	    window.open( link, \'popup\', \'width=900,height=700,resizable=yes,scrollbars=yes\' );
	}
        function fue( gene ) {
	    link = \'$PLOT_LINK2\' + gene + '.png';
	    window.open( link, \'popup\', \'width=900,height=700,resizable=yes,scrollbars=yes\' );
	}
	var newWin;
	function popupWindow( term, geneList, genes ) {
	    if( newWin != null ) {
			newWin.close();
		}
        newWin = window.open(\'http://yahoo.com\',\'popup\',\'width=900,height=700,resizable=yes\'); 
    }
    </script></head><body>

<table>
<tr><td valign=top>
<table border = 1>
<tr><th>FOLD_CHANGE_COLOR_KEY</th></tr>
<tr><td style=\"background-color:#6464f0\">Fold change &lt 0.25 </td></tr>
<tr><td style=\"background-color:#9292f0\"> 0.25 =&lt Fold change &lt 0.5</td></tr>
<tr><td style=\"background-color:#c1c1f0\"> 0.5 =&lt Fold change &lt 0.8</td></tr>

<tr><td style=\"background-color:#f0c1c1\"> 1.2 =&lt Fold change &lt 2</td></tr>
<tr><td style=\"background-color:#f09292\"> 2 =&lt Fold change &lt 4</td></tr>
<tr><td style=\"background-color:#f06464\"> 4 &lt Fold change</td></tr>

<tr><td style=\"background-color:f0f0f0\"> Other </td></tr>
</table>
</td>
<td valign=top>
<table border = 1>
<tr><th>pvalue_COLOR_KEY</th></tr>
<tr><td style=\"background-color:#FF0000\">p-value &lt 1.0E-4 </td></tr>
<tr><td style=\"background-color:#FF9700\"> 1.0E-4 =&lt p-value &lt 1.0E-3</td></tr>

<tr><td style=\"background-color:#FFCA00\"> 1.0E-3 =&lt p-value &lt 0.01</td></tr>
<tr><td style=\"background-color:#FFfd00\"> 0.01 =&lt p-value &lt 0.1</td></tr>
<tr><td style=\"background-color:#C0C0C0\"> other</td></tr>
</table>
</td>
<tr></table>

<h1><a href=\"$link\"><center>$TITLE : $pathway_name</center></a></h1>
<script type=\"text/javascript\" src=\"http://artemis.na.pg.com/sorttable.js\"></script>
<table border= 1 class=\"sortable\" id=\"93275610\">
";
    return $header;
}


sub time_string{
    my @Weekdays = ('Sunday', 'Monday', 'Tuesday', 'Wednesday', 'Thursday', 'Friday', 'Saturday');
    my @Months = ('January', 'February', 'March', 'April', 'May', 'June', 'July', 'August', 'September', 'October', 'November', 'December');
    my @Now = localtime(time());
    my $Month = $Months[$Now[4]];
    my $Weekday = $Weekdays[$Now[6]];
    my $Hour = $Now[2];
    my $AMPM;
    if ($Hour > 12) {
	$Hour = $Hour - 12;
	$AMPM = "PM";
    } else {
	$Hour = 12 if $Hour == 0;
	$AMPM ="AM";
    }
    my $Minute = $Now[1];
    $Minute = "0$Minute" if $Minute < 10;
    my $Year = $Now[5]+1900;
    my $time_string="$Weekday, $Month $Now[3], $Year $Hour:$Minute $AMPM";
    return $time_string;
}

sub fill_stat_filter{
    my ($ph,$fh,$file, $filter)=@_;
    my %id;
    open (X, "$filter")||die "could not open $filter\n";
    while (my $line = <X>){
	chomp $line;
	my @tmp=split /\t/, $line;
	$id{$tmp[0]}=1;
    }
    close X;
    open (X, "$file")||die "could not open $file\n";
    while (my $line = <X>){
	chomp $line;
	my @tmp=split /\t/, $line;	
	if (defined $id{$tmp[0]}){
	    if (@tmp<2){next;}
	    $$ph{$tmp[0]}=$tmp[1];
	    $$fh{$tmp[0]}=$tmp[2];
	}
    }
    close X;
    return;
}

sub fill_stat_log_filter{
    my ($ph,$fh,$file, $filter)=@_;
    my %id;
    open (X, "$filter")||die "could not open $filter\n";
    while (my $line = <X>){
	chomp $line;
	my @tmp=split /\t/, $line;
	$id{$tmp[0]}=1;
    }
    close X;
    open (X, "$file")||die "could not open $file\n";
    my $tmpx=<X>;
    while (my $line = <X>){
	chomp $line;
	my @tmp=split /\t/, $line;
	if (defined $id{$tmp[0]}){
	    if (@tmp<2){next;}
	    $$ph{$tmp[0]}=$tmp[1];
	    $$fh{$tmp[0]}=2**($tmp[2]);
	}
    }
    close X;
    return;
}
