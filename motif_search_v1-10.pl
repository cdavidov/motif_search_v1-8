#!/usr/bin/perl -w
#
# Chen Davidovich,  chen.davidovich@monash.edu

use strict;
use warnings;
use File::Basename;
use 5.010;
use Getopt::Long qw(GetOptions);

# If a fasta file for reference genome is not available, use combine_fasta_genome.pl to download the required build and to combine the chromosome fasta files to one fasta file. 

# Motif to identify:
#	(XnNm)k (
# The actual search for boundaries: (XnNm)(k-1)Xn, 
# to make sure that the first and last registers are for base of type X.
# E.g, if X_F=G, X_R=C, k=2, n=3, m=1-9:
# CGTGACGTCAGGGAAGCCCGCAGGGCCCG
# motif (F):[------------>]
# motif (R):      [<---------]
# motif (R):       [<--------]
  
my $sVarName = "";
my $sOutPath = ".";
my $sInputFastaFile = "";
 



my $n_min=1;#1,2,3,4,5
my $n_max='';# For not limiting the maximum write ''
my $m_min=1;
my $m_max=1;# For not limiting the maximum write '' --> WARRING: in case of m_max this is usually makes no sense, and in case of a large genome this will take almost forever and will create huge output file.
my $k_min=1;#1,2,3,4,5
my $k_max=1;# For not limiting the maximum write ''
my $sX_fwd="";#G,C,A,T
my $sX_rev="";#C,G,T,A
my $sN_fwd="[A|C|T|G|a|c|t|g|N|n]"; # e.g. [A|C|T|G] is any base
my $sN_rev="[T|G|A|C|t|g|a|c|N|n]"; # e.g. [A|C|T|G] is any base
#my $sN="[A|T|G]"; # e.g. [A|C|T|G] is any base
my $fMaxGC = 1.0;
my $fMinGC = 0;
my $sGC_regex = "GgCc";
my $bCheck_GC = 0;
my $bWrite_GC = 0;
my $bRegex_X=0;
my $bPrintOutInfoWhileRunning = 0;
my $bHelp = 0;
my $sModuleLoadCom = "";
my $sBedtoolsRunCom = "bedtools merge";
my $sBedtoolsOptions = "";
my $bMerge=0;
my $sManPage = 
'
Description:
	
	This script search within a genome (fasta file) for 
	a sequence motif that defined as follows:
		(XnNm)k 
	Where X and N are whildcard character(s) (perl format) 
	that can be defined seperatly for the forwared (X_fwd,N_fwd)
	or reverse (X_rev, N_rev) strands. n, m and k are integers.
	
	The script look for the bounderies of the motif where
	the actual boundaries are as follows: (XnNm)(k-1)Xn, 
	to make sure that the first and last registers 
	are for base/sequence of type X.
	
	Output bed file will create in the output directory.
	The bed file will include chromosome, start and end for each motif, 
	the sequence and the strand (-/+) from where it was identifed. 
	
	E.g, if X_fwd=[G|g], X_rev=[C|c], k=2, n=3, m=1-9 
	and N_fwd/N_rev=[A|C|T|G|a|c|t|g|N|n]:
	
	AACGTGACGTCAGGGAAGCCCGCAGGGCCCG
	motif 1 (F):[------------>]
	motif 2 (R):      [<---------]
	motif 3 (R):       [<--------]
	
	A seperate script, combine_fasta_genome.pl, can 
	be used to download the required build and to 
	combine the individual chromosome fasta files to 
	one fasta file.
	
Usage: 

	perl motif_search_vx-x.pl [OPTIONS]

Mandatory arguments:

	--fasta-in|-f
		Input genome in FASTA format.
		
	--X-fwd
		Motif, as will be searched on the + strand. Use 
		text in the format of perl whildcard string.
		e.g. if -X-fwd can be either G or g:
		--X-fwd \'[G|g]\'
				
	--X-rev
		Motif, as will be searched on the - strand. Use 
		text in the format of perl whildcard string.
		e.g. if -X-rev can be either C or c:
		--X-fwd \'[C|c]\'
		
Optional arguments:
	
	--var-name|-v = string
		A string that will be added to
		output file names (should be linux 
		compatible). 
		[Default: ""]
	
	--out-dir|o=s = string
		Directory to where output files will
		be written (no need to add the final \\).
		[Default: "."]
	
	--N-fwd
		Linker between repeated motifs, as will be searched 
		on the + strand. Use text in the format of perl 
		whildcard string.
		e.g. if -N-fwd can be any base:
		--N-fwd \'[A|C|T|G|a|c|t|g|N|n]\'
		[Default: \'[A|C|T|G|a|c|t|g|N|n]\']
		
	--N-rev
		Linker between repeated motifs, as will be searched 
		on the - strand. Use text in the format of perl 
		whildcard string.
		e.g. if -N-rev can be any base:
		--N-fwd \'[A|C|T|G|a|c|t|g|N|n]\'
		[Default: \'[T|G|A|C|t|g|a|c|N|n]\']
	
	--n-min
		Min n
		[Default: 1]
		
	--n-max
		Max n. Integer to specify a vaue or \'\' for infinity.
		[Default: \'\', i.e. infinity].
	
	--m-min
		Min m
		[Default: 1]
		
	--m-max
		Max m. For not limiting the maximum write \'\'.
		Yet, beware that unlimited --m-max is usually makes 
		no sense, and in case of a large genome will take 
		long time and will create huge output file.
		[Default: 1].
	
	--k-min
		Min k
		[Default: 1]
		
	--k-max
		Max k. Integer to specify a vaue or \'\' for infinity.
		[Default: \'\', i.e. infinity].	
		
	--max-GC
		Max GC content per motif (0.0 to 1.0)
		[Default: 1.0]
		
	--min-GC
		Min GC content per motif (0.0 to 1.0)	
		[Default: 1.0]
	
	--write-GC
		Add GC content into the output bed file
		(Unless --max-GC and min-GC set to other 
		value than the default, this will increase the run time)
		[Default: false]
	
	--regex-GC
		Characters that should be considered as
		either G or C for GC content calculation.
		[Default: GgCc]
		
	--regex-X
		If this flag is used, the strings passed through
		--X-fwd and --X-rev will be used as Perl regular expression
		arguments, without bringint into account m and k (n can still
		be used to determine how many time the motif will repeat).
		
	--merge
		use bedtools merge to merge overlapped motifs in the 
		output bed file (will create a second bed file, 
		named *.merge.bed)
		
	--module-load|-l=string
		System command used to load the module for 
		bedtools run (do not include if irrelevant).
	
	--bedtools-merge-com = string
		Command used to run bedtools merge. 
		[Default: \'bedtools merge\']
	
	--bedtools-merge-opt = string
		A string with optional arguments that will 
		be passed to bedtools closest. 
		[Default: \'\']
		
	--print
		print out info while running.
	
	--help
		Print this man page.
		
Examples:

	perl ./scripts/motif_search_vx-x.pl \
	-v test1 \
	-o ./MyData/Projects/GnNmSearch/2016m05d29-test \
	-f ./MyData/Projects/GnNmSearch/MotifSearch/ReferenceGenome/mm9/mm9.fa \
	--n-min 4 \
	--n-max \'\' \
	--m-min 1 \
	--m-max 4 \
	--k-min 4 \
	--k-max \'\' \
	--X-fwd \'[G|g]\' \
	--X-rev \'[C|c]\' \
	--print

';

GetOptions(
	"fasta-in|f=s" => \$sInputFastaFile, 
	"var-name|v=s" => \$sVarName,
	"out-dir|o=s" => \$sOutPath,
	"n-min=n" => \$n_min,
	"n-max=s" => \$n_max,# string, to allow for '' in case of no input.
	"m-min=n" => \$m_min,
	"m-max=s" => \$m_max,# string, to allow for '' in case of no input.
	"k-min=n" => \$k_min,
	"k-max=s" => \$k_max,# string, to allow for '' in case of no input.
	"X-fwd=s" => \$sX_fwd,
	"X-rev=s" => \$sX_rev,
	"N-fwd=s" => \$sN_fwd,
	"N-rev=s" => \$sN_rev,
	"max-GC=f" => \$fMaxGC,
	"min-GC=f" => \$fMinGC,
	"write-GC" => \$bWrite_GC,
	"regex-GC=s" => \$sGC_regex,
	"regex-X" => \$bRegex_X,
	"merge" => \$bMerge,
	"module-load|l=s" => \$sModuleLoadCom,
	"bedtools-merge-com=s" => \$sBedtoolsRunCom,
	"bedtools-merge-opt=s" => \$sBedtoolsOptions,
	"print" => \$bPrintOutInfoWhileRunning,
	"help" => \$bHelp
	) or die "Invalid input argument in $0\n";


if ($bHelp){
	print($sManPage."\n");
	exit 0;
}	
	
if ($sInputFastaFile eq ""){
	print ("ERROR: --fasta-in|-f is a mandatory argument. See manual (--help).\nTerminating.\n");
	exit 0;
}

if ($sX_fwd eq ""){
	print ("ERROR: --X-fwd is a mandatory argument. See manual (--help).\nTerminating.\n");
	exit 0;
}

if ($sX_rev eq ""){
	print ("ERROR: --X-rev is a mandatory argument. See manual (--help).\nTerminating.\n");
	exit 0;
}
if (($fMaxGC<1) or ($fMinGC>0)){ #user would like to check GC
	$bCheck_GC = 1;
} 


#--------------------------------------------------------

if (-e $sOutPath){
	# directory is already present, do nothgin.
 }else{
	unless (mkdir $sOutPath){ die "Unable to create $sOutPath\n. Either there is no writing permition to the parental directory or it might not be exist.";}
 }	


if ($bPrintOutInfoWhileRunning){
	my $datestring = localtime();
	print "Staet time: $datestring\n";
	print "\$n_min=$n_min; ";
	print "\$n_max=$n_max; ";# For not limiting the maximum write ''
	print "\$m_min=$m_min; ";
	print "\$m_max=$m_max; ";# For not limiting the maximum write '' --> WARRING: in case of m_max this is usually makes no sense, and in case of a large genome this will take almost forever and will create huge output file.
	print "\$k_min=$k_min; ";#
	print "\$k_max=$k_max; ";# For not limiting the maximum write ''
	print "\$sX_fwd=$sX_fwd; ";
	print "\$sX_rev=$sX_rev; ";
	print "\$sN_fwd=$sN_fwd; "; # e.g. [A|C|T|G] is any base
	print "\$sN_rev=$sN_rev; "; # e.g. [A|C|T|G] is any base
	print "\$fMinGC=$fMinGC; ";
	print "\$fMaxGC=$fMaxGC; ";
	print "\$sGC_regex=$sGC_regex; ";
	print "\$bRegex_X=$bRegex_X; ";
	print "\$bMerge=$bMerge; ";
	print "\$sModuleLoadCom=$sModuleLoadCom; ";
	print "\$sBedtoolsRunCom=$sBedtoolsRunCom; ";
	print "\$sBedtoolsOptions=$sBedtoolsOptions; ";
	print "\n";
	}

my ($name,$path);
($name,$path) = fileparse($sInputFastaFile);

my $sX_fwd_clean = $sX_fwd;
$sX_fwd_clean =~ s/[^[:alnum:]_-]//g;

my $sN_fwd_clean = $sN_fwd;
$sN_fwd_clean =~ s/[^[:alnum:]_-]//g;
 
my $sOutFile = $sOutPath."/motifs_".$sVarName."_fasta=$name.bed";


my @ahFastaData = @{&getFasta($sInputFastaFile)};
my (@aPos);
my ($raStart_fwd,$raEnd_fwd);
my ($raStart_rev,$raEnd_rev);
my (@aStart_fwd);
my (@aEnd_fwd);
my (@aStart_rev);
my (@aEnd_rev);


if ($bPrintOutInfoWhileRunning){ print "Open $sOutFile for writing.\n";}
open (my $fhBedOut, '>', $sOutFile) or die "can't open file $sOutFile for writing\n";

my $k_min_minus_one = $k_min - 1;
my $k_max_minus_one;
if ($k_max eq ''){
		$k_max_minus_one = '';
	}else{
		$k_max_minus_one = $k_max - 1;
	}

	my $sMotif_fwd;
	my $sMotif_rev;
	if($bRegex_X){
		$sMotif_fwd = "((".$sX_fwd."){".$n_min.",".$n_max."})";
		$sMotif_rev = "((".$sX_rev."){".$n_min.",".$n_max."})";
	}else{
		$sMotif_fwd = "(".$sX_fwd."{".$n_min.",".$n_max."}".$sN_fwd."{".$m_min.",".$m_max."}){".$k_min_minus_one.",".$k_max_minus_one."}".$sX_fwd."{".$n_min.",".$n_max."}";
		$sMotif_rev = "(".$sX_rev."{".$n_min.",".$n_max."}".$sN_rev."{".$m_min.",".$m_max."}){".$k_min_minus_one.",".$k_max_minus_one."}".$sX_rev."{".$n_min.",".$n_max."}";


		}


		
my $sMotif;
my $sStrand;
my $sIdentifiedMotif;
my $j_fwd;
my $j_rev;


my ($fGC_content,$nNumberOf_GC_bases,$nStringLength);

for (my $i=1;$i<scalar(@ahFastaData);$i++){ #iterate over all chromosomes
	
	
 			($raStart_fwd,$raEnd_fwd) = match_all_positions($sMotif_fwd,$ahFastaData[$i]->{'seq'});
			@aStart_fwd = @{$raStart_fwd};
			@aEnd_fwd = @{$raEnd_fwd};
			
			($raStart_rev,$raEnd_rev) = match_all_positions($sMotif_rev,$ahFastaData[$i]->{'seq'});
			@aStart_rev = @{$raStart_rev};
			@aEnd_rev = @{$raEnd_rev};
			
			$j_fwd=0;
			$j_rev=0;
			
			
			for (my $j=0;$j<(scalar(@aEnd_fwd)+scalar(@aEnd_rev));$j++){ # go over all end and start position of motifs  at teh fwd and rev strands and write them out by the order of their chromosomal position.
				
				if ((($j_fwd<scalar(@aEnd_fwd)) and ($j_rev==scalar(@aEnd_rev))) or (($j_fwd<scalar(@aEnd_fwd)) and ($j_rev<scalar(@aEnd_rev)) and (((scalar(@aEnd_rev)==0) and (scalar(@aEnd_fwd)>0)) or ($aStart_fwd[$j_fwd] <= $aStart_rev[$j_rev])))) {# use the fwd
					$sIdentifiedMotif = substr $ahFastaData[$i]->{'seq'}, ($aStart_fwd[$j_fwd] - 1) , $aEnd_fwd[$j_fwd] - $aStart_fwd[$j_fwd] + 1;
						
						
						$fGC_content=0;
						if (($bCheck_GC) or ($bWrite_GC)){
						($fGC_content,$nNumberOf_GC_bases,$nStringLength) = check_GC($sIdentifiedMotif,$sGC_regex);
						# print "mot: $sIdentifiedMotif; GC=$nNumberOf_GC_bases; fGC=$fGC_content\n";
						}
						
						#######
						#print $fhBedOut $ahFastaData[$i]->{'label'}."\$bCheck_GC=$bCheck_GC \$fGC_content=$fGC_content \$fMinGC=$fMinGC \$fMaxGC=$fMaxGC\t+\n";
						#######
						
						if (($bCheck_GC==0) or (($bCheck_GC) and ($fGC_content>=$fMinGC) and ($fGC_content<=$fMaxGC))){
							print $fhBedOut $ahFastaData[$i]->{'label'}."\t$aStart_fwd[$j_fwd]\t$aEnd_fwd[$j_fwd]\t$sIdentifiedMotif\t$fGC_content\t+\n";
							#$j_fwd++;
							#print "(\$fGC_content,\$nNumberOf_GC_bases,\$nStringLength)=($fGC_content,$nNumberOf_GC_bases,$nStringLength)\n"
						}
						$j_fwd++;
					
				}
				
				if ((($j_fwd==scalar(@aEnd_fwd)) and ($j_rev<scalar(@aEnd_rev))) or (($j_rev<scalar(@aEnd_rev)) and ($j_fwd<scalar(@aEnd_fwd)) and (((scalar(@aEnd_fwd)==0) and (scalar(@aEnd_rev)>0)) or ($aStart_fwd[$j_fwd] > $aStart_rev[$j_rev])))){
					$sIdentifiedMotif = substr $ahFastaData[$i]->{'seq'}, ($aStart_rev[$j_rev] - 1) , $aEnd_rev[$j_rev] - $aStart_rev[$j_rev] + 1;
					
					#############
					    $fGC_content=0;
						if (($bCheck_GC) or ($bWrite_GC)){
						($fGC_content,$nNumberOf_GC_bases,$nStringLength) = check_GC($sIdentifiedMotif,$sGC_regex);
						 #print "mot: $sIdentifiedMotif; GC=$nNumberOf_GC_bases; fGC=$fGC_content\n";
						}
						
						#############
						if (($bCheck_GC==0) or (($bCheck_GC) and ($fGC_content>=$fMinGC) and ($fGC_content<=$fMaxGC))){
							print $fhBedOut $ahFastaData[$i]->{'label'}."\t$aStart_rev[$j_rev]\t$aEnd_rev[$j_rev]\t$sIdentifiedMotif\t$fGC_content\t-\n";
						
							#$j_rev++;
						}
						$j_rev++;
					
				}
				
				
				
			}#j
	 
}#i

	 
close $fhBedOut;

# merge the bed file, if required:
my $sCom;
if ($bMerge){ #user would like to merge the bed file:
	
	$sCom = "$sModuleLoadCom; $sBedtoolsRunCom $sBedtoolsOptions -i $sOutFile > $sOutFile.merge.4col.bed";
	if ($bPrintOutInfoWhileRunning){print ($sCom."\n");}#print
	system($sCom);#run
	
	#$sCom = 'awk \'{print $1,"\t",$2,"\t",$3,"\t",".","\t","0","\t",$4}\''." $sOutFile.merge.4col.bed > $sOutFile.merge.bed";
	$sCom = 'awk \'{print $1"\t"$2"\t"$3"\t.\t0\t"$4}\''." $sOutFile.merge.4col.bed > $sOutFile.merge.bed";
	
	if ($bPrintOutInfoWhileRunning){print ($sCom."\n");}#print
	system($sCom);#run

}


if ($bPrintOutInfoWhileRunning){
	my $datestring = localtime();
	print "End time: $datestring\n";
}



########################## subs ##############
sub getFasta(){
# Get sequences and descriptors (labels) from a fasta file with multiple sequences.

#E.g.:
# my @ahFastaData = @{&getFasta("file_name.fa")};
# print $ahFastaData[1]->{'seq'}."\n"; #Print the first sequence.

	my $sInputFastaFile = $_[0];
	#my $fhLogOut = $_[1];
	my @ahFastaData;

	my $iSequence=0;
	if ($bPrintOutInfoWhileRunning){ print "open file $sInputFastaFile for reading.\n";}
	open (my $fhFastaIn,'<',$sInputFastaFile) or die "can't open file $sInputFastaFile for reading\n";
	my $sSequenceDescriptor;
	my $bStartingNewSequence = 0;
	my $sSequence;
	my $nNewLineLength = 0;
	my $nSeqLength = 0;
	while ( my $line = <$fhFastaIn> ) {
		chomp($line);
		$nNewLineLength = length($line);
		if ($bStartingNewSequence==1){ #this line starts a new sequence (i.e., next line after a line started with ">")
			
			
			
			$ahFastaData[$iSequence]->{'seq'} = $line;
			$bStartingNewSequence = 0;
		}elsif ($line =~ m/^>(.*)/){ #This is a descriptor line, starts with ">"
			
			$iSequence++;
			
			$ahFastaData[$iSequence]->{'label'} = $1;
			$bStartingNewSequence = 1;
			if ($bPrintOutInfoWhileRunning){ print "Read data for field: $ahFastaData[$iSequence]->{'label'}\n";}
		}else{ # this is a sequence line, which is not the first sequence line in this sequence.
			$nSeqLength = length($ahFastaData[$iSequence]->{'seq'});
			
			substr($ahFastaData[$iSequence]->{'seq'}, $nSeqLength, 0, $line); # FAST: at 50 M string: 5 sec/10%-increment less than 0.1 sec fast 
			
			}
	}
	
	close $fhFastaIn;
	return \@ahFastaData;
}#sub

sub match_all_positions {
    
	# Get a text string and a regular expression string as text, and return two arrys with start and stop positions of all matched sub strings.
	# e.g."
	# ($raStart,$raEnd) = match_all_positions("GGG",$ahFastaData[$i]->{'seq'});
	# @aStart = @{$raStart};
	# @aEnd = @{$raEnd};
	
    my ($regex, $string) = @_;
    my @reg_start;
	my @reg_end;
    #while ($string =~ /$regex/g) {
	
	while ($string =~ /(?=($regex))/g) {
        # The ?= is to ensure that overlapping motifs will be captured.
				
		push @reg_start, ($-[0] + 1);
		push @reg_end, $+[1];
    }
    return (\@reg_start,\@reg_end);
}

# e.g.:
#my ($fGC_content,$nNumberOf_GC_bases,$nStringLength) = check_GC($sequence,$sGC_regex);
sub check_GC{

	my ($string,$sGC_regex) = @_;
	my $nNumberOf_GC_bases;
	eval('$nNumberOf_GC_bases = $string =~ tr/'.$sGC_regex.'//;');
	my $fGC_content = $nNumberOf_GC_bases/length($string);
	#print "$string (\$fGC_content,\$nNumberOf_GC_bases,\$nNumberOf_GC_bases)=($fGC_content,$nNumberOf_GC_bases,$nNumberOf_GC_bases)\n";
	return ($fGC_content,$nNumberOf_GC_bases,length($string))
					

}

