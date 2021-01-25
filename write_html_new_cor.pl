#! /usr/bin/perl -w
##########################################################################
#### Author: Ping HU
#### I need to rewrite the program and make things much simpler, 
#### one theme just one gene page
#### also each theme showed whether it is up or down-regulation 
#### like the combine data
###Get hash of
########1. stat
########2. ann syn and long
########3. present file
########4. ontorank term and edge
###1. step: run ontorank, sig, up, down
###2. combine data into minP, allP, upP, downP, allN, upN, downN, allGeneList
###3. Display data
###4. Display all data at the top link to title
###5. can consider seperate BP, MF and CL in future to add a row so to seperate some process
##########################################################################
use strict;
use Getopt::Long;
use List::Util qw(first max maxstr min minstr reduce shuffle sum);


my $PLOT_LINK2="http://mvic-biotech.na.pg.com/projects/GoldenEagle/GSS2527/Bioinformatics/boxplot/";

my $PLOT_LINK="http://mvic-biotech.na.pg.com/projects/GoldenEagle/GSS2527/Bioinformatics/boxplot/";


my ($opt_l);
GetOptions("l=s" => \$opt_l);
my $usage = "$0 player_list annfile P_CUT list_file filter_file
   -l chromosome location file
   \n";

die $usage unless @ARGV == 7;
my ($player_list,$annfile, $P_CUT, $file, $long_name_file, $acr_file, $filter) = @ARGV;

my $PLAYF= $player_list;
if ($player_list =~  /([^\/]+)\s*$/ ){
    $PLAYF=$1;
}
my $PLAYN=get_player_N($player_list);

my $SUBTITLE= $file;
if ($file =~  /([^\/]+)\s*$/ ){
    $SUBTITLE=$1;
}
if ($SUBTITLE =~/(\S+)\.txt/){
    $SUBTITLE=$1;
}
my %LONG_NAME;
my %acr;
fill_hash($long_name_file,\%LONG_NAME);
fill_hash($acr_file,\%acr);


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

my %PPP;
my %FFF;
#####Cut data into list to run the autorank2 #########
my ($allsiggene, $allsigN)=cut_list_filter($file, $SUBTITLE, \%PPP, \%FFF, $filter);

my $name_all=$SUBTITLE."_p05";
my $name_up=$SUBTITLE."_upp05";
my $name_dn=$SUBTITLE."_dnp05";

system ("/home/ping/pkg/ontorank_java/new_onto_rank/bin/runOntoRank2 -r $name_all /home/ping/db/affy/NA36/GO.edge  $annfile $player_list /home/ping/db/affy/NA36/GO.term  > $name_all.result\n");
system ("/home/ping/pkg/ontorank_java/new_onto_rank/bin/runOntoRank2 -r $name_up /home/ping/db/affy/NA36/GO.edge  $annfile $player_list  /home/ping/db/affy/NA36/GO.term > $name_up.result\n");
system ("/home/ping/pkg/ontorank_java/new_onto_rank/bin/runOntoRank2 -r $name_dn /home/ping/db/affy/NA36/GO.edge  $annfile $player_list /home/ping/db/affy/NA36/GO.term  > $name_dn.result\n");

########Combine result together#######################
my %data;
my %CN;
my %TotalGN;
my %Genelist;
open (A, "<$name_all.result")||die "could not open $name_all.result\n";
while (my $a=<A>){
    chomp $a;
    my @tmp=split /\t/, $a;
    if(scalar @tmp >=5){
	my $PP=$tmp[0]."\t".$tmp[5];
	$data{$PP}{"all"}=$tmp[4];
	$CN{$PP}{"all"}=$tmp[1];
	$TotalGN{$PP}=$tmp[2];
	$Genelist{$PP}=$tmp[6];
    }
}
close A;

open (A, "<$name_up.result")||die "could not open $name_up.result\n";
while (my $a=<A>){
    chomp $a;
    my @tmp=split /\t/, $a;
    if(scalar @tmp >=5){
	 my $PP=$tmp[0]."\t".$tmp[5];
	 $data{$PP}{"up"}=$tmp[4];
	 $CN{$PP}{"up"}=$tmp[1];
    }
    
}
close A;

open (A, "<$name_dn.result")||die "could not open $name_dn.result\n";
while (my $a=<A>){
    chomp $a;
    my @tmp=split /\t/, $a;
    if(scalar @tmp >=5){
	 my $PP=$tmp[0]."\t".$tmp[5];
	 $data{$PP}{"dn"}=$tmp[4];
	 $CN{$PP}{"dn"}=$tmp[1];
    }
    
}
close A;

########Now need to write the result to a html page
my $outfile="$SUBTITLE.html";
system("mkdir ./$SUBTITLE.html_files");

my $sub_item_number=0;
open (B, ">$outfile")or die "could not open $outfile\n ";
my $MAIN_HEADER= get_main_header($SUBTITLE);
print B $MAIN_HEADER;
$sub_item_number++;
my $sub_item_file_name="./$SUBTITLE.html_files/ITEM".$sub_item_number.".html";
write_gene_html ($allsiggene, $sub_item_file_name, "TotalPresentProbes: $PLAYN ---TotalSignificantProbe: $allsigN--- $SUBTITLE", \%PPP, \%FFF, \%LONG_NAME, \%acr);
print B "<h2><a href=\"$sub_item_file_name\" target=\"_geneWindow\" );>Significant Gene $allsigN in Total $PLAYN present probes </a></h2>\n";
print B "<table border=1 class=\"sortable\" id=\"5429378998752\" >\n<tr><td>Theme</td><td>trend</td><td>minP</td><td>P.all</td><td>SigNum.all</td><td>P.up</td><td>SigNum.up</td><td>P.down</td><td>SigNum.down</td><td>TotalProbeNumber</td><tr>\n";
my %outline;
foreach my $i(sort keys %data){
    my $genestring=$Genelist{$i};
    my $totalGN=$TotalGN{$i};
    
    for($i){s/\s*null\s*//gi;}
    my $trend="--";
    my $color="black";
    my $pall=1;

    my $nall=0;
    if(defined $data{$i}{"all"}){
	$pall=$data{$i}{"all"};
	$nall=$CN{$i}{"all"};
    }

    my $pup=1;
    my $nup=0;
     if(defined $data{$i}{"up"}){
	$pup=$data{$i}{"up"};$nup=$CN{$i}{"up"};
    }
    my $pdn=1;my $ndn=0; 
    if(defined $data{$i}{"dn"}){
	$pdn=$data{$i}{"dn"};$ndn=$CN{$i}{"dn"};
    }
 
    my $minP=min($pall, $pup, $pdn);
    if(($minP <=0.05)&&($nall>=3) &&(defined $genestring)&&(!($genestring =~ /^\s*$/))){
	$trend="sig";
	$color="green";
	if($pdn < $pup){$trend="down"; $color="blue";}
	if($pdn > $pup){$trend="up";$color="red";}

	$sub_item_number++;
	my $sub_item_file_name="./$SUBTITLE.html_files/ITEM".$sub_item_number.".html";
	if($nup==0){
	    $nup=$nall-$ndn;
	}elsif($ndn==0){
	    $ndn=$nall-$nup;
	}
	write_gene_html ($genestring, $sub_item_file_name, "$i ---trend: $trend ---ThemeP= $minP ---Totalprobes: $totalGN ---TotalSignificantProbe: $nall (up: $nup down: $ndn )--- $SUBTITLE", \%PPP, \%FFF, \%LONG_NAME, \%acr);
	$outline{$minP}{$sub_item_file_name}= "<tr style=\"color: $color\"><td>". $i. "</td><td><a href=\"$sub_item_file_name\" target=\"_geneWindow\" );>". $trend."</a></td><td>".$minP."</td><td>" .$pall. "</td><td>".$nall. "</td><td>".  $pup. "</td><td>". $nup. "</td><td>". $pdn. "</td><td>". $ndn. "</td><td>".$totalGN. "</td></tr>\n";
    }
}
foreach my $i (sort {$a<=>$b} keys %outline){
    foreach my $j (keys %{$outline{$i}}){
	print B $outline{$i}{$j};
    }
}
my $tail=get_end();
print B  "\n</table>", $tail;
close B;


sub fill_hash{
    my ($file, $HF)=@_;
    open (X, "<$file")or die "could not open $file\n ";
    while (my $line=<X>){
        chomp $line;
        my @tmp=split /\t/, $line;
        if (defined $tmp[1]){
            $$HF{$tmp[0]}=$tmp[1];
        }else{
            $$HF{$tmp[0]}="#";
        }
    }
    close X;
    return;
}


sub get_span_P{
    my ($p, $P_CUT)=@_;
    my $span_class="style=\"background-color:C0C0C0\"";
    my $score=0;
    if ($p <=0.001*$P_CUT){
	$span_class="style=\"background-color:FF0000\"";
	$score=4;
    }elsif($p<=0.01*$P_CUT){
	$span_class="style=\"background-color:FF9700\"";
	$score=3;
    }elsif($p<=0.1*$P_CUT){
	$span_class="style=\"background-color:FFca00\"";
	$score=2;
    }elsif($p<=$P_CUT){
	$span_class="style=\"background-color:FFfd00\"";
	$score=1;
    }

    return ($span_class, $score);
}




sub get_color_code{
    my ($P_CUT)=@_;
    my $content="<table border=\"1\" align=\"CENTER\">\n"
    . "<tr><td align=\"CENTER\"><h2>Heat Map Color Legend</h2></td></tr>\n<tr>\n<td style=\"background-color:C0C0C0\">\n"
    .$P_CUT
    ." < p</td>\n</tr>\n<tr>\n<td style=\"background-color:FFfd00\">"
    .0.1*$P_CUT
    . " < p <= "
    .$P_CUT
    ."</td>\n</tr>\n<tr>\n<td style=\"background-color:FFca00\">"
    .0.01*$P_CUT
    ." < p <= "
    .0.1*$P_CUT
    ."</td>\n</tr>\n<tr><td style=\"background-color:FF9700\">"
    .0.001*$P_CUT
    . " < p <= "
    .0.01*$P_CUT
    ."</td>\n</tr>\n<tr>\n<td style=\"background-color:FF0000\">p <= "
    .0.001*$P_CUT
    ."</td>\n</tr>\n</table>\n";
    return $content;
}

sub get_end{
    my $tail="<h3>Report created on: ".time_string()."</h3>\n</body>\n</html>\n";
    return $tail;
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



sub write_gene_html{
    ######gene names will be #seperated###########
    my ($gene_names, $file_name, $tname, $HPPP, $HFFF, $HLONG_NAME,$Hacr) =@_;
    my @array =split /\#/, $gene_names;
    open (DDDD, ">$file_name")or die "could not open $file_name\n ";

    my $header =get_sub_header($tname);
    print DDDD $header;
    print DDDD "\n\n<table border=1 class=\"sortable\" id=\"5429378998758\" >\n<tr>\n";
    print DDDD "<td><font color=\"#000099\">Gene</font></td>\n";
    print DDDD "<td><font color=\"#000099\">Acronym</font></td>\n";
    print DDDD "<td><font color=\"#330033\">P</font></td>\n";
    print DDDD "<td><font color=\"#996600\">Fold</font></td>\n";
    print DDDD "<td><font color=\"#006600\">Name</font></td>\n";
    if ($opt_l){
	print DDDD "<td><font color=\"#006600\">Location</font></td>\n";	
    }
    print DDDD "</tr>\n";
    my %order;my $other="";
    for (my $i =0; $i < scalar @array; $i++){
	my $name=$array[$i];
	my $P="N/A";
	my $FC="N/A";
	my $long="N/A";
	my $aaa="N/A";
	if (defined $$HPPP{$name}){$P=$$HPPP{$name};}
	if (defined $$HFFF{$name}){$FC=$$HFFF{$name};}
	if (defined $$HLONG_NAME{$name}){$long=$$HLONG_NAME{$name};}
	if (defined $$Hacr{$name}){$aaa=$$Hacr{$name};}

	my $span_class="style=\"background-color:f0f0f0\"";
	my $span_p="style=\"background-color:C0C0C0\"";
	if (($FC ne "N/A")&&($FC ne "")&&($FC ne "-")){
	    if ($FC >0.5){
		$span_class="style=\"background-color:f06464\""; 
	    }elsif($FC>0.3){
		$span_class="style=\"background-color:f09292\"";
	    }elsif($FC>0){
		$span_class="style=\"background-color:f0c1c1\"";
	    }elsif($FC==0){
		$span_class="style=\"background-color:f0f0f0\"";
	    }elsif($FC>-0.3){
		$span_class="style=\"background-color:c1c1f0\"";
	    }elsif($FC>-0.5){
		$span_class="style=\"background-color:9292f0\"";
	    }else{
		$span_class="style=\"background-color:6464f0\"";
	    }
	}
	if (($P ne "N/A")&&($P ne "")&&($P ne "-")){
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
       
	my $line="<tr>\n<td><a href=\"javascript:void(0)\" onclick=\"gar(\'$aaa\');\">$name</a> </td>\n"
	    ."<td><a href=\"javascript:void(0)\" onclick=\"gsr(\'$aaa\');\">$aaa</a> </td>\n"
	    ."<td $span_p><a href=\"javascript:void(0)\" onclick=\"fue(\'$aaa\');\">$P</a></td>\n"."<td $span_class>$FC</td>\n"
	    ."<td><a href=\"javascript:void(0)\" onclick=\"fue(\'$aaa\');\"><font color=\"#006600\">$long</font></a></td>\n"
	    ."</tr>\n";
	
	if ($opt_l){
	    my $loc="N/A";
	    if (defined $location{$name}){
		$loc= $location{$name}
	    } 
	    $line="<tr>\n<td><a href=\"javascript:void(0)\" onclick=\"gar(\'$aaa\');\">$name</a> </td>\n"
	    ."<td><a href=\"javascript:void(0)\" onclick=\"gsr(\'$aaa\');\">$aaa</a> </td>\n"
	    ."<td $span_p><a href=\"javascript:void(0)\" onclick=\"fue(\'$aaa\');\">$P</a></td>\n"."<td $span_class>$FC</td>\n"
	    ."<td><a href=\"javascript:void(0)\" onclick=\"fue(\'$aaa\');\"><font color=\"#006600\">$long</font></a></td>\n"
	    ."<td><font color=\"#006600\">$loc</font></td>\n"
	    ."</tr>\n";
	}
	if (($FC ne "N/A")&&($FC ne "")&&($FC ne "-")){
	    if (! defined $order{$FC}){$order{$FC}=$line;}
	    else{
		$order{$FC}.=$line;
	    }
	}else{
	    $other.=$line;
	} 
    }
    foreach my $i (sort {$a<=>$b} keys  %order){
	print DDDD $order{$i};
    }
    print DDDD $other;
    print DDDD "</table>\n"; 
    print DDDD time_string(), "\n";
    print DDDD "</BODY></HTML>\n";
    return;
}

sub get_main_header{
    my ($TITLE)=@_;
    my $MAIN_HEADER="<html>\n<head>\n
    <title>Clusterer Heat Map</title>
	<script language=\"JavaScript\">
	function oldgsr( gene ) {
	    link = \'http://en.wikipedia.org/wiki/\'+ gene;
		window.open( link, \'popup\', \'width=900,height=700,resizable=yes,scrollbars=yes\' );
	}
	function gsr( gene ) {
	    link = \'http://www.genecards.org/cgi-bin/carddisp.pl?gene=\'+ gene;
		window.open( link, \'popup\', \'width=900,height=700,resizable=yes,scrollbars=yes\' );
	}
	function gar( gene ) {
	    link = \'$PLOT_LINK\' + gene + '.Derm.png';
	    window.open( link, \'popup\', \'width=900,height=700,resizable=yes,scrollbars=yes\' );
	}
        function fue( gene ) {
	    link = \'$PLOT_LINK2\' + gene + '.Epi.png';
	    window.open( link, \'popup\', \'width=900,height=700,resizable=yes,scrollbars=yes\' );
	}

	var newWin;
	function popupWindow( term, geneList, genes ) {
	    if( newWin != null ) {
			newWin.close();
		}
        newWin = window.open(\'http://yahoo.com\',\'popup\',\'width=900,height=700,resizable=yes\'); 
    }

	</script>
        <script type=\"text/javascript\" src=\"http://artemis.na.pg.com/sorttable.js\"></script>

  </head>

  <body>
  <h1><center>$TITLE</center></h1>
 
";
    return $MAIN_HEADER;
}

sub get_sub_header{
    my ($name)=@_;
        my $header="<html>
	<head>

	<script language=\"JavaScript\">
	
	function oldgsr( gene ) {
	    link = \'http://en.wikipedia.org/wiki/\'+ gene;
		window.open( link, \'popup\', \'width=900,height=700,resizable=yes,scrollbars=yes\' );
	}
	function gsr( gene ) {
	    link = \'http://www.genecards.org/cgi-bin/carddisp.pl?gene=\'+ gene;
		window.open( link, \'popup\', \'width=900,height=700,resizable=yes,scrollbars=yes\' );
	}


	function gar( gene ) {
	    link = \'$PLOT_LINK\' + gene + '.Derm.png';
	    window.open( link, \'popup\', \'width=900,height=700,resizable=yes,scrollbars=yes\' );
	}
        function fue( gene ) {
	    link = \'$PLOT_LINK2\' + gene + '.Epi.png';
	    window.open( link, \'popup\', \'width=900,height=700,resizable=yes,scrollbars=yes\' );
	}



	var newWin;
	function popupWindow( term, geneList, genes ) {
	    if( newWin != null ) {
			newWin.close();
		}
        newWin = window.open(\'http://yahoo.com\',\'popup\',\'width=900,height=700,resizable=yes\'); 
    }

	</script>
<script type=\"text/javascript\" src=\"http://artemis.na.pg.com/sorttable.js\"></script>
  </head>

  <body>
  
  <h1><center>$name<center></h1>

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

";
    return $header;
}

sub cut_list_filter{
    ##generate the sublist with $SUBTITLE list file, of p from 0.05, 0.1, 0.01 
    my ($file, $SUBTITLE, $P_hash, $F_hash, $filterfile)=@_;
    my %id;
    my $siggene="";
    my $sigN=0;
    open (XXXX, "$filterfile")||die "could not open $filterfile\n";
    while (my $line = <XXXX>){
	chomp $line;
	my @tmp=split /\t/, $line;
	$id{$tmp[0]}=1;
    }
    close XXXX;
    open (XXXX, "<$file")or die "could not open $file\n ";
    my $tmp=<XXXX>;
    while (my $line=<XXXX>){
	chomp $line;
	my @tmp=split /\t/, $line;
	if (defined $id{$tmp[0]}){
	    if (defined $tmp[1]){
		$$P_hash{$tmp[0]}=$tmp[1];
	    }
	    if (defined $tmp[2]){
		$$F_hash{$tmp[0]}=$tmp[2];
	    }
	}
    }
    close XXXX;
      
    my $name4=$SUBTITLE."_p05";
    my $name4U=$SUBTITLE."_upp05";
    my $name4D=$SUBTITLE."_dnp05";
    my $name=$SUBTITLE;
    
   
    open (A4, ">$name4")or die "could not open $name4\n "; 
    open (A4U, ">$name4U")or die "could not open $name4U\n ";
    open (A4D, ">$name4D")or die "could not open $name4D\n ";
   
    
    foreach my $i (keys %PPP){


	if ($PPP{$i} <=0.05){
	    $sigN++;
	    if($siggene eq ""){
		$siggene=$i;
	    }else{
		$siggene =$siggene."#". $i;
	    }
	    print A4 $i, "\n";
	    if ($FFF{$i}>0){
		print A4U $i, "\n";
	    }elsif($FFF{$i}<0){
		print A4D $i, "\n";
	    }
	}

    }
 
   
    close A4;
   
    close A4U;
   
    close A4D;
    
    open (A, ">$name")or die "could not open $name\n ";
  
  
    print A $name4, "\n";
   
    print A $name4U, "\n";
    
    print A $name4D, "\n";
   
    close A;
    return($siggene, $sigN);
    
}
    
sub get_player_N{
    my ($player_list)=@_;
    open (XX, "<$player_list")or die "could not open $player_list\n ";
    my $PLAYN=0;
    while (my $line =<XX>){
	chomp $line;
	if (!($line =~ /^\s*$/)){
	    $PLAYN++;
	}
    }
    close XX;
    return $PLAYN;
}
sub get_gene_string{
    my ($player_list)=@_;
    my $result="";
    open (XXX, "<$player_list")or die "could not open $player_list\n ";
    while (my $line =<XXX>){
	chomp $line;
	if (!($line =~ /^\s*$/)){
	    my @tmp=split /\t/, $line;
	    if ($result eq ""){
		$result.=$tmp[0];
	    }else{
		$result.="#".$tmp[0];
	    }
	}
    }
    close XXX;
    return $result;
}


sub get_long_name{
    my ($long_name_file, $Hacr, $HLONG_NAME) =@_;
open (X1, "<$long_name_file")or die "could not open $long_name_file\n ";
while (my $line=<X1>){
    my @tmp=split /\t/, $line;
    if (defined $tmp[2]){
	$$HLONG_NAME{$tmp[0]}=$tmp[2];
    }
    if (defined $tmp[1]){
	$$Hacr{$tmp[0]}=$tmp[1];
    }        
}
close X1;
return;
}

