#!/usr/bin/perl
use strict;
use Digest::SHA1 qw(sha1_hex);  # sha1_hex 関数をインポートしておく

&main;

sub main{
    my %database;
    %database=&read_sequence;
    &show_result(\%database);
}

sub read_sequence{
    my $line;
    my $name="";
    my $sequence="";
    my $sequence_f="";
    my $number=0;
    my %database;
    my %duplicate;
    my $identical=0;
    
    while($line=<stdin>){
	$line=~s/[\r\n]//g;
	if($line=~/^>/){
	    if($name ne ""){
		if(exists($database{$sequence})){
		    $duplicate{$database{$sequence}[0]}.="$name\n";
		    $identical++;
		}else{
		    $database{$sequence}[0]=$name;
		    $database{$sequence}[1]=$sequence_f;
		    $database{$sequence}[2]=$number;
		}
	    }
	    $name=$line;
	    $sequence="";
	    $sequence_f="";
	    $number++;
	}else{
	    $sequence.=$line;
	    $sequence_f.=$line."\n";
	}
    }
    if($name ne ""){
	if(exists($database{$sequence})){
	    $duplicate{$database{$sequence}[0]}.="$name\n";
	    $identical++;
	}else{
	    $database{$sequence}[0]=$name;
	    $database{$sequence}[1]=$sequence_f;
	    $database{$sequence}[2]=$number;
	}
    }
    &show_duplicates(\%duplicate);
    print "# $identical sequences are identical to the others\n";
    return %database;
}

sub show_result{
    my %database=%{$_[0]};
    my $key;
    foreach $key (sort {$database{$a}[2] <=> $database{$b}[2]} keys %database){
#	print "#$database{$key}[2]\n";
	print $database{$key}[0]."\n";
	print $database{$key}[1];
#	print "#\n";
    }
}

sub show_duplicates{
    my %duplicates=%{$_[0]};
    my $key;
    foreach $key (keys %duplicates){
	print "# $key is identical with ...\n";
	$duplicates{$key}=~s/^/# /g;
	$duplicates{$key}=~s/\n/\n# /g;
	print "$duplicates{$key}";
	print "\n";
    }
}
