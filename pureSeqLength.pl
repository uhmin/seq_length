#!/usr/bin/perl
use strict;
our $seq_length="/home/uhmin/bin/seq_length.pl";
&main;

sub main{
    my %options;
    my %region;
    my @result;
    my $preN;

    %options=&interface;
    if(!exists($options{-i})){
	$options{-i}=&createtmpfile;
	$options{delete}=$options{-i};
    }
    while(1){
	@result=&calc_region(\%options);
	if($result[0] == $preN){
	    print stderr "$preN sequences were left.\n";
	    &output(\%options);
	    last;
	}elsif($result[0] < 2){
	    print stderr "No sequence was left with this parameter.\n";
	    last;
	}
	$preN=$result[0]*1;
    }
    if(exists($options{delete})){
	unlink($options{delete});
    }
    return 0;
}

sub interface{
    my %options;
    my $i;
    $options{-w}="1";
    for($i=0; $i<@ARGV; $i+=2){
	if($ARGV[$i] eq "-i" || $ARGV[$i] eq "-d" 
	   || $ARGV[$i] eq "-r" || $ARGV[$i] eq "-w"){
	    $options{$ARGV[$i]}=$ARGV[$i+1];
	}else{
	    print stderr "Unknown option: $ARGV[$i]\n";
	    &help();
	}
    }
    return %options;
}

sub help{
    my $file=__FILE__;
    print stderr << "EOF";
$file: Set distribution of sequence length within specified standard deviation.
  -i infile
  -d standard deviation; more than 3 is recommended.
  -r region; If you want to specify length limit use this option.
  -w width
EOF
;
    exit;
}

sub calc_region{ #(\%options);
    my $command="";
    my $infile=$_[0]->{-i};
    my $line;
    my @result;
    my $from;
    my $to;

    if(exists($_[0]->{-r})){
	$command="$seq_length -i $infile -r $_[0]->{-r} | ";
    }else{
	$command="cat $infile | ";
    }
    $command.="$seq_length";
#    print "$command\n";
    open FIN, "$command | " or die("Could not exec the command:\n $command");
    while($line=<FIN>){
#	print $line;
	if($line=~/^number of sequence = /){
	    ($result[0]=$line)=~s/.* = //;
	}
	if($line=~/^average\s+= /){
	    ($result[1]=$line)=~s/.* = //;
	}
	if($line=~/^variance = /){
	    ($result[2]=$line)=~s/.* = //;
	}
    }
    close FIN;

    printf(stderr "%d sequences, average: %.2f , variance: %.2f  "
	   , $result[0], $result[1], $result[2]);
    $from=int($result[1]-$result[2]*$_[0]->{-d}+.5);
    if($from<1){
	$from=1;
    }
    $to  =int($result[1]+$result[2]*$_[0]->{-d}+.5);
    $_[0]->{-r}="$from-$to";
    print stderr "Set region: $_[0]->{-r}\n";
    return @result;
}

sub output{#(\%options);
    my $command;
    my $result;
    $command="$seq_length -i $_[0]->{-i} -r $_[0]->{-r}";
    $result=`$command | $seq_length -w $_[0]->{-w}`;
    print stderr $result;
    print `$command`;
    return 0;
}

sub createtmpfile{
    my $basename="tmp";
    my $i=0;
    my $filename;
    my $line;
    $filename=$basename;
    while( -r $filename ){
	$filename="${basename}${i}";
	$i++;
    }
    open FOUT, ">$filename" 
	or die("Could not create temporaly file: $filename\n");
    while($line=<stdin>){
	print FOUT $line;
    }
    close FOUT;
    return $filename;
}
