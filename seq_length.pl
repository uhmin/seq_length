#! /usr/bin/perl
use strict;
our $WIDTH=80;
our $prename="";
our $line;

&main;


sub main{
    my %result;
    my %option;
    my %histogram;
    my %statistics;
    $result{"fileend"}=1;
    %option=&interface;
    &heading;
    undef %histogram;
    if($option{"mode"} eq "histogram"){
        %statistics=&make_histogram(\%histogram, \%option);
        &display_statistic(\%statistics);
        &display_histogram(\%histogram, \%statistics, \%option);
    }elsif($option{"mode"} eq "sequence"){
        while($result{"fileend"}){
            %result=&read_single_fasta;
            if($option{"from"}<=$result{"length"}
               && ($result{"length"}<=$option{"to"}
                   || $option{"to"}==0)){
                print $result{"seqname"};
                if($option{"add_length"} eq "T"){
                    print " length=".$result{"length"};
                }
                print "\n";
                print $result{"sequence"};
            }
        }
    }
    exit;
}

sub heading {
    while($line=<FIN>){
        #print"* $line";
        if($line=~/^>/){
            $prename=$line;
            chomp($prename);
            last;
        }
    }
}

#%statistics=&make_histogram(\%histogram, \%option);
sub make_histogram{
# main_rutin
    my %statistics;
    my $length;
    my $keylength;
    my %result;
    my $windowsize=$_[1]->{"windowsize"};

    $result{"fileend"}=1;
    $statistics{"sum"}=0;
    $statistics{"sum2"}=0;
    $statistics{"number"}=0;

        %result=&read_single_fasta;
        $keylength=int($result{"length"}/$windowsize)*$windowsize;
        $_[0]->{"$keylength"}++;
        $statistics{"sum"}+=$result{"length"};
        $statistics{"sum2"}+=$result{"length"}**2;
        $statistics{"number"}++;
    $statistics{"min"}=$result{"length"};
    $statistics{"max"}=$result{"length"};

    while($result{"fileend"}){
        %result=&read_single_fasta;
        $keylength=int($result{"length"}/$windowsize)*$windowsize;
        $_[0]->{"$keylength"}++;
        $statistics{"sum"}+=$result{"length"};
        $statistics{"sum2"}+=$result{"length"}**2;
        $statistics{"number"}++;
        if($statistics{"max"} < $result{"length"}){
            $statistics{"max"} = $result{"length"};
        }
        if($statistics{"min"} > $result{"length"}){
            $statistics{"min"} = $result{"length"};
        }
    }
    return %statistics;
}

sub read_single_fasta{
    my %result;
    my $length;
    my $sequence;
    my $seqname;
    $sequence="";
    $seqname=$prename;
    chomp($line);
    if(!$line=~/^>/){$length=length($line);}
    else{$length=0;}
    while($line=<FIN>){
        if($line=~/^>/){
            $result{"fileend"}=1;
            $prename=$line;
            chomp($prename);
            last;
        }
        else{
            $sequence.=$line;
            chomp($line);
            $length+=length($line);
        }
        $result{"fileend"}=0;
    }
    $result{"seqname"}=$seqname;
    $result{"length"}=$length;
    $result{"sequence"}=$sequence;
    return %result;
}


sub display_statistic{
    my $sum=$_[0]->{"sum"};
    my $sum2=$_[0]->{"sum2"};
    my $number=$_[0]->{"number"};
    my $min=$_[0]->{"min"}*1;
    my $max=$_[0]->{"max"}*1;
    print FOUT "number of sequence = ".$number."\n";
    if($number>1){
        print FOUT "length: $min-$max\n";
        print FOUT "average  = ".($sum/$number)."\n";
        print FOUT "variance = ".sqrt(($sum2-$sum*$sum/$number)/($number-1))."\n";
    }
}

#&display_histogram(\%histogram, \%statistics, \%options);
sub display_histogram{
    my $i;
    my @list;
    my $number;
    my @max;
    my $key;
    my %histogram=%{$_[0]};
    my $space;
    my $result;
    my $cumlative=0;
    my @keta;
    my $format;

    @list=sort{$a <=> $b} keys %histogram;

    $max[2]=0;
    foreach $key(keys %histogram){
        if($max[2]<$histogram{$key}){$max[2]=$histogram{$key}}
    }

    $number=scalar(@list);
    $max[0]=$list[0];
    $max[1]=$list[$number-1];
    if($number==0 || $max[1]<=0 || $max[2]<=0){
        print FOUT "Input file seem to be wrong.\n";
        return 0;
    }
    $keta[0]=log($max[1])/log(10)+1;
    $keta[1]=log($_[1]->{"number"})/log(10)+1;
    $keta[2]=log($max[2])/log(10)+1;
    $format=sprintf("%%%dd-%%-%dd %%%dd %%%dd (%%5.1f) "
                    , $keta[0], $keta[0], $keta[2], $keta[1]);
    $WIDTH-=length($format);

    foreach $key (sort {$a <=> $b} keys %histogram){
        $cumlative+=$histogram{$key};
        $result=sprintf($format,
                        , $key, $key+$_[2]->{"windowsize"}-1
                        , $histogram{$key}, $cumlative, $cumlative*100/$_[1]->{"number"});
        if($_[2]->{"windowsize"}<2){
            $result=~s/.*-//;
        }
        print FOUT $result;
        for($i=0;$i<$WIDTH*$histogram{$key}/$max[2];$i++){
            print FOUT "*";
        }
        print FOUT "\n";
    }
}


sub interface(){
    my $i;
    my %option;
    my $from;
    my $to;

    $option{"mode"}="histogram";
    if(@ARGV % 2 !=0){
        &help;
    }
    *FIN=*STDIN;
    *FOUT=*STDOUT;
    $option{"windowsize"}=1;
    foreach ($i=0;$i<@ARGV;$i+=2){
        if($ARGV[$i] eq "-w"){
            $option{"windowsize"}=$ARGV[$i+1];
        }elsif($ARGV[$i] eq "-i"){
            open FIN, "$ARGV[$i+1]" or &list_open_err($ARGV[$i+1]);
        }elsif($ARGV[$i] eq "-o"){
            open FOUT, ">$ARGV[$i+1]" or &list_open_err($ARGV[$i+1]);
        }elsif($ARGV[$i] eq "-r"){
            $option{"mode"}="sequence";
            ($option{"from"},$option{"to"})=split(/-/,$ARGV[$i+1]);
        }elsif($ARGV[$i] eq "-l"){
            ($option{"add_length"}=$ARGV[$i+1])=~s/(\w+)/\U$1\E/g;
        }else {
            print $ARGV[$i]."\n";
            &help();
        }
    }
    return %option;
}


sub help{
    print stderr __FILE__."\n";
    print stderr "Output histogram of sequence length from multifasta.\n";
    print stderr "\nOPTIONS:\n";
    print stderr " -i [infile name(%s)] stdin if default.\n";
    print stderr " -o [outfile name(%s)] stdout id default.\n";
    print stderr " -w [windowsize of the histogram(%d)]. One if default.\n";
    print stderr " -r [from(%d)]-[to(%d)] extract sequence length of [from] to [to]. if to==0 lengsh limit is infinit.\n";
    print stderr " -l [T/F] add sequence length after sequence name\n";
    exit;
}

sub list_open_err{
    print stderr "Could not open file ".$_[0]."\n";
    exit;
}
