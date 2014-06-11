#!/usr/bin/perl
use strict;

&main;

sub main{
    my $line;
    my $locus;
    my $definition;
    my $accession;
    my $gi;
    my $organism;
    my $country="NoCountry";
    my $genotype="NoGenotype";
    my $date="NoDate";
    my $isolate="";
    my $sequence="";
    my $tmp;
    
    my $sw;
    for($sw=0; $line=<stdin>;){
	chomp($line);
	if($sw==1){
	    if($line=~/^\/\//){
		print ">${genotype}|${country}|${date}|$gi|$accession|$locus $isolate $definition $organism\n";
		print "$sequence\n";
		$sw=0;
		$locus="";
		$genotype="NoType";
		$country="NoCountry";
		$date="NoDate";
		$isolate="";
		$organism="";
		$sequence="";
		$accession="";
		$gi="";
	    }else{
		$line=~s/^\s*\d+\s+//;
		$line=~s/\s+//g;
		$sequence.=$line."\n";
	    }
	}elsif($line=~/^LOCUS\s/){
	    ($locus=$line)=~s/LOCUS\s+(\S*)\s.*/$1/;
#	    print "locus: $locus\n";
	}elsif($line=~/^DEFINITION/){
	    ($definition=$line)=~s/^DEFINITION\s+//;
	}elsif($line=~/^ACCESSION/){
	    ($accession=$line)=~s/ACCESSION\s+//;
	    $accession=~s/\s.*//;
	}elsif($line=~/^VERSION/ && $line=~/GI:/){
	    ($gi=$line)=~s/.*GI://;
	    $gi=~s/\s.*//;
	}elsif($line=~/\/organism=/){
	    ($organism=$line)=~s/.*\/organism="(.*)"/$1/;
#	    print "organism: $organism\n";
	}elsif($line=~/\/country=/){
	    ($country=$line)=~s/.*\/country="(.*)"?/$1/;
	    $country=~s/:.*//;
	    $country=~s/"//g;
	    $country=~s/\s+/_/g;
#	    print "country: $country\n";
	}elsif($line=~/\/note="/){
	    while($line!~/\/note=".*"/){
		$tmp=<stdin>;
		$tmp=~s/^\s+/ /;
		chomp($tmp);
		$line.=$tmp;
	    }
	    if($line=~/genotype/i){
		($genotype=$line)=~s/.*genotype(.*)"/$1/i;
		$genotype=~s/^\W*//;
		$genotype=~s/\s+/_/g;
#	    print "genotype: $genotype\n";
	    }
	}elsif($line=~/\/collection_date=/i){
	    ($date=$line)=~s/\/collection_date="(.*)"/$1/i;
	    $date=~s/.*-//;
	    $date=~s/\s+//g;
#	    print "collection date: $date\n";
	}elsif($line=~/\/isolate=/i){
	    ($isolate=$line)=~s/\/isolate="(.*)"/$1/i;
	    $isolate=~s/^\s+//;
	    $isolate=~s/\s+/_/g;
#	    print "collection date: $date\n";
	}elsif($line=~/^ORIGIN/){
	    $sw=1;
	    $sequence="";
	}elsif($line=~/genotype\W/ && $genotype eq "NoType"){
	    ($genotype=$line)=~s/.*genotype//;
	    $genotype=~s/^\W+(\S+)\s.*/$1/;
	    $genotype.="(?)";
	}
    }
    if($locus ne ""){
	print ">${genotype}|${country}|${date}|$locus $organism\n";
	print "$sequence\n";
    }
}
