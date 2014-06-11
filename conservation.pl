#!/usr/bin/perl
use strict;

&main;

sub main{
    my $seqnum;
    my @analysisResult;
    my @data;
    my %in;
    %in=&interface;
    @analysisResult=&read_sequence(\%in);
    @data=&analyze(@analysisResult);
    if($in{m} eq "s"){
	&display_sequence($analysisResult[0], \@data, 1);
    }elsif($in{m} eq "x"){
	&display_tabular($analysisResult[0], \@data);
    }else{
	&display_sequence($analysisResult[0], \@data, 0);
    }
}

sub interface{
    my $i; my $key;
    my $argc=scalar(@ARGV);
    my %in;
    if($argc % 2 ==1){ &help; }
    foreach($i=0; $i<$argc; $i++){
        if($ARGV[$i]=~/-/){
            ($key=$ARGV[$i])=~s/^-//;
            $in{$key}=$ARGV[$i+1];
            if($key ne "i" && $key ne "m"){ # -i infile -m mode
                print stderr "unknown option: $key\n";
                &help;
            }
        }
    }
    return %in;
}

sub read_sequence{
    my %in=%{$_[0]};
    my $line;
    my  $i=-1;
    my $j;
    my @seqname;
    my @sequence;
    my @length;
    my @result;
    $main::maxlength=0;
    if(exists($in{i})){
	open FIN, $in{i} or die;
    }else{
	*FIN=*STDIN;
    }
    while($line=<FIN>){
	if($line=~/^>/){ #seqname;
	    $i++;
	    $seqname[$i]=$line;
	    $sequence[$i]="";
	}else{
	    $line=~s/[\r|\n|\s]//g;
	    $sequence[$i].=$line;
	}
    }
    close FIN;
    for($j=0; $j<=$i; $j++){
	$length[$j]=length($sequence[$j]);
	if($main::maxlength<$length[$j]){
	    $main::maxlength=$length[$j];
	}
    }
    $result[0]=$i;
    $result[1]=\@sequence;
    $result[2]=\@seqname;
    $result[3]=\@length;
    return @result;
}

sub analyze{
    my $seqnum=$_[0];
    my @sequence=@{$_[1]};
    my @length=@{$_[3]};
    my $i; my $j;
    my $sw=1;
    my @data;
    my $base;
    for($j=0; $sw ;$j++){
	$sw=0;
	for($i=0; $i<=$seqnum; $i++){
	    if($i==0){
		undef(%{$data[$j]});
	    }
	    if($j<$length[$i]){
		$base=substr($sequence[$i],$j,1);
		$sw=1;
	    }else{
		$base="-";
	    }
	    $data[$j]{$base}+=1;
	    $main::baselist{$base}=1;
	}
    }
    return @data;
}

sub display_sequence{
    my $seqnum=$_[0];
    my @data=@{$_[1]};
    my $i;
    my $length=@data;
    my $sw;
    my $result;
    my $all_base;
    my $base;

    print ">Input ".($seqnum+1)." sequences\n";
    for($i=0; $i<$main::maxlength; $i++){
#    for($i=0; $i<$length; $i++){
	$sw=0;
	$result="";
	$all_base="";
	foreach $base (sort {$data[$i]{$b} <=> $data[$i]{$a}} keys %{$data[$i]}){
	    if($sw!=0){
		$result.=", ";
	    }
	    $result.=$base.":".$data[$i]{$base};
	    $sw++;
	    $all_base.=$base;
	    if($_[2] == 0){
		last;
	    }
	}
	if($sw==1){
	    print $all_base;
#	    if($all_base ne "-"){
#		print $all_base;
#	    }
	}elsif($all_base ne "-"){
	    print "($result)";
	}
    }
    print "\n";
}

sub display_tabular{
    my $seqnum=$_[0];
    my @data=@{$_[1]};
    my $i;
    my $base;
    my $length=@data;
    my $sw;
    my @baselist_local;
    my @max;
    my $result;
    my $rate;
    my $uncertainty;

    @baselist_local=keys(%main::baselist);

    print ">Input ".($seqnum+1)." sequences\n";
    print "#";
    foreach $base (@baselist_local){
	print "\t".$base;
    }
    print "\n";

    for($i=0; $i<$main::maxlength; $i++){
#    for($i=0; $i<$length; $i++){
	$sw=0;
	$result=$i+1;
#	print $i;
	$max[0]=0;$max[1]="-";
	$uncertainty=0;
	foreach $base (@baselist_local){
	    if(exists($data[$i]{$base})){
		$rate=$data[$i]{$base}/($seqnum+1);
		$uncertainty-=$rate*(log($rate)/log(2));;
		$result.="\t".$data[$i]{$base};
		$sw+=$data[$i]{$base};
#		print "\t".$data[$i]{$base};
		if($max[0]<$data[$i]{$base}){
		    $max[0]=$data[$i]{$base};
		    $max[1]=$base;
		}
	    }else{
		$result.="\t0";
#		print "\t0";
	    }
	}
	if($sw>0){
	    $result.=sprintf("\t%s\t%5.1f%%"
			     , $max[1], int($max[0]/($seqnum+1)*1000)/10);
	    $result.="\t$uncertainty";
	    $result.="\n";
	    print $result;
	}
#	printf("\t%s(%5.1f%%)", $max[1], int($max[0]/($seqnum+1)*1000)/10);
#	print "\n";
    }
    print "\n";
}

sub help{
    print stderr __FILE__."\n";
    print stderr "Usage: \n";
    print stderr "-i infile\n";
    print stderr "-m output mode [s/x]. s: fasta like format, x: tabular format.\n";
    exit;
}
