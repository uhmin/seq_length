#!/usr/bin/perl
use strict;

&main;

sub main{
    my @region;
    my %options;
    my @result;
    %options=&interface();
    @result=&readData(\%options);
}

sub interface{
    my %options;
    my $i;
    my $key;
    my $value;
    $options{"-i"}="F";
    $options{"-a"}="F";
    for($i=0; $i<@ARGV; $i++){
	$key=$ARGV[$i];
	$value=$ARGV[$i+1];
	if($key eq "-i" || $key eq "-a" 
	   || $key eq "-c"){
	    $options{$key}=$value;
	    $i++;
	}elsif($key eq "--inverse"){
	    $options{"-i"}="T";
	}elsif($key eq "--anysubject"){
	    $options{"-a"}="T";
	}else{
	    print stderr "Unknown option: $key\n";
	    &help;
	}
    }
    if($options{"-i"} ne "T" && $options{"-i"} ne "F"){
	print stderr "-i option should be T or F\n";
	&help;
    }
    if($options{"-a"} ne "T" && $options{"-a"} ne "F"){
	print stderr "-a option should be T or F\n";
	&help;
    }
    return %options;
}

sub help{
    my $file=__FILE__;
    $file=~s/.*\///;
    print stderr << "EOF";

-c: [integer] Cuotff cover ratio. 
    $file will discard subjects having 
    less than this cover ratio.

-a: [T/F] If you specify this option, 
    $file will ignore subject sequence name.
    This means, $file will calculate cover region 
    of any subject sequence at once.

-i: [T/F] If you specify this option, $file will calculate 
    the cover region of the subject sequence.

--anysubject: Same meaning as -a T
--inverse:    Same meaning as -i T

  This program calculates coverd region of query sequence.
$file < blast_tabular_file > output

  If coverd region of subject regopm, do this.
$file --inverse < blast_tabular_file > output
$file -i T      < blast_tabular_file > output

  If you want to calculate cover region of such 
   as alot of small sequences to the query sequence, do this.
$file --anysubject < blast_tabular_file > output
$file -a T         < blast_tabular_file > output

EOF
;
    exit;
}

sub readData{ #(%options);
    my $options=$_[0];
    my $line; 
    my $i=0;
    my @data; 
    my @tmpdata;
    my @region;
    my $blindDB=$options->{"-a"};

    while($line=<stdin>){
	if($line!~/\w/){ next; }
	chomp($line);
	@tmpdata=split(/\t/,$line);
	if($options->{"-i"} eq "T"){
#	    $,="\t";
#	    print @tmpdata;
#	    print "\n";
	    @tmpdata=&inverse(@tmpdata);
#	    print @tmpdata;
#	    print "\n";
	}

	if($blindDB eq "T"){
	    $tmpdata[1]="all";
	}
	if($i>0 && 
	   ($tmpdata[0] ne $data[$i-1][0] || $tmpdata[1] ne $data[$i-1][1])
	    ){
#	    print "$data[$i][0] v.s. $data[$i-1][0]\n";
#	    print "$data[$i][1] v.s. $data[$i-1][1]\n\n";
	    @region=&process(@data);
	    &showResult(\@region, \@data, $options);
	    $i=0;
	    @data=();
#	    @{$data[$i]}=split(/\t/,$line);
#	    if($blindDB==1){
#		$data[$i][1]="all";
#	    }
	}
	@{$data[$i]}=@tmpdata;
	$i++;
    }
    @region=&process(@data);
    &showResult(\@region, \@data, $options);
}

sub inverse{ #(@tmpdata);
    my @data=@_;
    my $swap;
    $swap=$data[0];
    $data[0]=$data[1];
    $data[1]=$swap;

    $swap=$data[6];
    $data[6]=$data[8];
    $data[8]=$swap;

    $swap=$data[7];
    $data[7]=$data[9];
    $data[9]=$swap;

    $swap=$data[12];
    $data[12]=$data[13];
    $data[13]=$swap;


    if($data[7] < $data[6]){
	$swap=$data[6];
	$data[6]=$data[7];
	$data[7]=$swap;

	$swap=$data[8];
	$data[8]=$data[9];;
	$data[9]=$swap;
    }
    return @data;
}

sub process{
    my $i; my $j;
    my @region=();
    my $numRegion=1;
    my $expand;

    my @tmp=@_;
    my @data=sort{${$a}[6] <=> ${$b}[6] 
		      || ${$a}[7] <=> ${$b}[7]} (@tmp);
    my $dataLength=scalar(@data);
    $region[0][0]=$data[0][6];
    $region[0][1]=$data[0][7];

    for($i=1; $i<$dataLength; $i++){
#	print "#".$data[$i][6]."\t".$data[$i][7]."\n";
	for($j=0, $expand=0; $j<$numRegion; $j++){
	    if($data[$i][6]<=$region[$j][1]+1){
		$region[$j][1]=max($data[$i][7], $region[$j][1]);
		$expand=1;
	    }
	}
	if(!$expand){
#	    print "# Expand region.\n";
	    $region[$numRegion][0]=$data[$i][6];
	    $region[$numRegion][1]=$data[$i][7];
	    $numRegion++;
	}
    }
    return @region;    
}

sub max{
    my $max=$_[0];
    my $tmp;
    foreach $tmp (@_){
	if($tmp>$max){
	    $max=$tmp;
	}
    }
    return $max;
}

sub showResult{
    my @region=@{$_[0]};
    my @data=@{$_[1]};
    my $options=$_[2];
    my @result;
    my $i; my $j; my $regionNum=scalar(@region);
    my $hitLength; my $sumLength;
    my $preposition=0;
    for($i=0, $sumLength=0; $i<$regionNum; $i++){
	$result[$i][0]=$data[0][0];
	$result[$i][1]=$data[0][1];
	$hitLength=$region[$i][1]-$region[$i][0]+1;
	if($data[0][12]==0){
	    $result[$i][2]="-";
	}else{
	    $result[$i][2]=$hitLength*100/$data[0][12];
	}
	$result[$i][3]=$hitLength;;
	$result[$i][4]="-";
	$result[$i][5]="-";
	$result[$i][6]=$region[$i][0];
	$result[$i][7]=$region[$i][1];
	$result[$i][8]=$preposition;
	$result[$i][9]=$region[$i][0]-1;
	$result[$i][10]="-";
	$result[$i][11]="-";
	$result[$i][12]=$data[0][12];
	$result[$i][13]=$data[0][13];
	$result[$i][14]=$data[0][14];
	$result[$i][15]=$data[0][15];
	$result[$i][16]=$data[0][16];
	$result[$i][17]="";
	$sumLength+=$hitLength;
	$preposition=$region[$i][1]+1;
    }
    if($sumLength*100/$data[0][12] >= $options->{"-c"}){
	&printResult($regionNum, \@data, $sumLength, \@result);
    }
}

sub printResult{
    my $regionNum=$_[0];
    my $data=$_[1];
    my $sumLength=$_[2];
    my $result=$_[3];
    my $i;
    my $j;
    for($i=0; $i<$regionNum; $i++){
	for($j=0; $j<17; $j++){
	    if($j==11){
		if($data->[0][12]==0){
		    printf("--\t");
		}else{
		    printf("%.1f\t", $sumLength*100/$data->[0][12]);
		}
	    }elsif($j==2){
		printf("%.1f\t",$result->[$i][$j]);
	    }elsif($j==10){
		print $sumLength."\t";
	    }else{
		print $result->[$i][$j]."\t";
	    }
	}
#	print $sumLength."----\n";
	print "\n";
    }
    print "\n";
}
