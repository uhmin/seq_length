#!/usr/bin/perl
use strict;

#our $fpicup="/user2/uhmin/bin_64/fpicup";
our $fpicup="fpicup";
our $head=10;
our $tail=10;
our $templete="Q";
our $identity=30;

&main;

sub main{
    my %options;
    %options=&interface;
    &read_data(\%options);
}

sub interface{
    my $i;
    my %options;
    my $sw=0;
    if(scalar(@ARGV) % 2 == 1){
	&help;
    }

    $options{-H}=$head;
    $options{-T}=$tail;
    $options{-e}=$templete;
    $options{-HC}="T";
    $options{-TC}="T";
    $options{-I}=$identity;

    for($i=0; $i<scalar(@ARGV); $i+=2){
	if($ARGV[$i] eq "-f" || $ARGV[$i] eq "-t"
	   || $ARGV[$i] eq "-H" || $ARGV[$i] eq "-T"
	   || $ARGV[$i] eq "-HC" || $ARGV[$i] eq "-TC"
	   || $ARGV[$i] eq "-e" || $ARGV[$i] eq "-I"){
	    $options{$ARGV[$i]}=$ARGV[$i+1];
	}else{
	    print stderr "Unknown option: $ARGV[$i]\n";
	    $sw=1;
	}
    }
    if($options{-e} ne "Q" && $options{-e} ne "D"){
	print stderr "-e option must be Q or D\n";
	$sw=1;
    }
    if($sw==1){
	&help;
    }
    return %options;
}

sub help{
    my $filename=__FILE__;
    print stderr << "EOF";
usage:
$filename -f fasta file name -t blast result table ( -H head gap -T tail gap)
 -f: Fasta sequence file
 -t: Blast result tabular format
 -H: Permissive length of gap at the start side. Default is $head.
 -T: Permissive length of gap at the  stop side. Default is $tail.
 -HC, -TC: [T/F] Whether if elongate the missing head or tail, respectively.
 -e: [Q/D] templete is Query or Database sequence. Defauls is $templete.
 -I: identity cut off. Default is $identity.
EOF
    ;
    exit;
}

sub read_data{ #(\%options);
    my $tablefile;
    my $fasta;
    my $line;
    my @data;
    my $from;
    my $to;
    my $key;
    my $command;
    my %duplicate;
    my $result;
    my @dataname;
    $tablefile=$_[0]->{-t};
    $fasta=$_[0]->{-f};
    open FIN, $tablefile or die("Could not open $tablefile\n");
    while($line=<FIN>){
	@data=split(/\t/, $line);
	if(scalar(@data)<13){
	    next;
	}
	if($data[2] < $_[0]->{-I}){
	    print stderr "Identity is less than threshold: $data[0] $data[1] $data[2]\n";
	    next;
	}

	if($_[0]->{-e} eq "D"){
	    &exchange(\@data);
	}

	&check_database(\@data);
	($from, $to)=&calc_region(\@data, $_[0]->{-H}, $_[0]->{-T}, $_[0]->{-HC}, $_[0]->{-TC});
	$key=$data[1];

	if($_[0]->{-e} eq "D"){
	    &exchange(\@data);
	}
	if($from!=0 && $to!=0){
#	    print $line;
	    if(exists($duplicate{$key}{$from-$to})){
		print stderr "Skip $key $from-$to: it already exists\n";
	    }else{
		$command="$fpicup -i $fasta -k \"$key\" -r $from-$to";
#		print $command."\n";
		$result= `$command`;
		if($result eq ""){
		    @dataname=split(/\|/, $key);
		    $key=$dataname[1];
		    $command="$fpicup -i $fasta -k \"$key\" -r $from-$to";
		    $result=`$command`;
		    if($result ne ""){
			print stderr "Found it as key: $key\n";
		    }
		}
		print $result;
	    }
	    $duplicate{$key}{$from-$to}=1;
	}else{
	    print stderr "not enough match: $key ";
	    print stderr "Q:$data[6]-$data[7]($data[12]) ";
	    print stderr "D:$data[8]-$data[9]($data[13])\n";
	}
    }
    close FIN;
}

sub check_database{ #(\@data);
    my $data=$_[0];
    my $tmp;
    my $direction=1;
    my $qlength;
    my $dlength;
    if($data->[6] > $data->[7]){ # qstart must be smaller than qstop
	$tmp=$data->[6];
	$data->[6]=$data->[7];
	$data->[7]=$tmp;
	
	$tmp=$data->[8];
	$data->[8]=$data->[9];
	$data->[9]=$tmp;
    }

#    print "# $data->[6]-$data->[7]  $data->[8]-$data->[9]\n";
    if($data->[8] > $data->[9]){
#	print "  Found Reverse hit  $data->[8]-$data->[9] ==>>  ";
	$data->[8]=&inverse($data->[8], $data->[13]);
	$data->[9]=&inverse($data->[9], $data->[13]);
	$direction=-1;
#	print "$data->[8]-$data->[9]   $data->[13]\n";
    }

    $qlength=$data->[7]-$data->[6]+1;
    $dlength=$data->[9]-$data->[8]+1;
    if($dlength*2.5 <= $qlength && $qlength <= $dlength*3.5){
#	print stderr "This data is from blastx!\n";
	# in case of blastx
	$data->[8]=$data->[8]*3-2;
	$data->[9]*=3;
	$data->[13]*=3;
    }

    if($qlength*2.5 <= $dlength && $dlength <= $qlength*3.5){
#	print stderr "This data is from tblastn!\n";
	# in case of blastx
	$data->[6]=$data->[6]*3-2;
	$data->[7]*=3;
	$data->[12]*=3;
    }
    push(@{$data}, $direction);
}

sub calc_region{ #(\@data, $_[0]->{-H}, $_[0]->{-T}, $_[0]->{-HC}, $_[0]->{-TC});
    my $from;
    my $to;
    my $limitHead=$_[1];
    my $limitTail=$_[2];
    my $HC=$_[3];
    my $TC=$_[4];
    my $data=$_[0];
    my $correction;
    my $sw=1;
    my $direction;
    my $correction;

    #head correction
    $correction=$data->[6]-1;
    if($correction > $limitHead){
	$from=0;
	$sw=0;
    }else{
	if($HC eq "F"){
	    $correction=0;
	}
#	print "correction:  $correction\n";
	$from=$data->[8]-$correction;
	if($from < 1){
	    $from=1;
	}
    }

    #tail correction
    $correction=$data->[12]-$data->[7];
    if($correction > $limitTail){
	$to=0;
	$sw=0;
    }else{
	if($HC eq "F"){
	    $correction=0;
	}
	$to=$data->[9]+$correction;
	if($to > $data->[13]){
	    $to=$data->[13];
	}
    }
    if($sw==0){
	$from=0;
	$to=0;
    }else{
	$direction=pop(@{$data});
	if($direction==-1){
#	    print "    Invert result $from-$to  ==>>";
	    $from=&inverse($from, $data->[13]);
	    $to  =&inverse($to,   $data->[13]);
#	    print "  $from-$to $data->[12]  --------------\n";
	}
    }
    return ($from, $to);
}

sub inverse{
    my $result;
    $result=$_[1]-$_[0]+1;
    return $result;
}

sub exchange{ #(\@data);
    my $tmp;
    my $data=$_[0];
    
    # name
    $tmp=$data->[0];
    $data->[0]=$data->[1];
    $data->[1]=$tmp;

    # start
    $tmp=$data->[6];
    $data->[6]=$data->[8];
    $data->[8]=$tmp;

    # stop
    $tmp=$data->[7];
    $data->[7]=$data->[9];
    $data->[9]=$tmp;

    #length
    $tmp=$data->[12];
    $data->[12]=$data->[13];
    $data->[13]=$tmp;
}
