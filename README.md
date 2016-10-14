## motif_search_v1-10

# Description:
	
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
	
# Usage: 

	perl motif_search_vx-x.pl [OPTIONS]
	
# Mandatory arguments:

	--fasta-in|-f
		Input genome in FASTA format.
		
	--X-fwd
		Motif, as will be searched on the + strand. Use 
		text in the format of perl whildcard string.
		e.g. if -X-fwd can be either G or g:
		--X-fwd '[G|g]'
				
	--X-rev
		Motif, as will be searched on the - strand. Use 
		text in the format of perl whildcard string.
		e.g. if -X-rev can be either C or c:
		--X-fwd '[C|c]'

# Optional arguments:
	
	--var-name|-v = string
		A string that will be added to
		output file names (should be linux 
		compatible). 
		[Default: ""]
	
	--out-dir|o=s = string
		Directory to where output files will
		be written (no need to add the final \).
		[Default: "."]
	
	--N-fwd
		Linker between repeated motifs, as will be searched 
		on the + strand. Use text in the format of perl 
		whildcard string.
		e.g. if -N-fwd can be any base:
		--N-fwd '[A|C|T|G|a|c|t|g|N|n]'
		[Default: '[A|C|T|G|a|c|t|g|N|n]']
		
	--N-rev
		Linker between repeated motifs, as will be searched 
		on the - strand. Use text in the format of perl 
		whildcard string.
		e.g. if -N-rev can be any base:
		--N-fwd '[A|C|T|G|a|c|t|g|N|n]'
		[Default: '[T|G|A|C|t|g|a|c|N|n]']
	
	--n-min
		Min n
		[Default: 1]
		
	--n-max
		Max n. Integer to specify a vaue or '' for infinity.
		[Default: '', i.e. infinity].
	
	--m-min
		Min m
		[Default: 1]
		
	--m-max
		Max m. For not limiting the maximum write ''.
		Yet, beware that unlimited --m-max is usually makes 
		no sense, and in case of a large genome will take 
		long time and will create huge output file.
		[Default: 1].

	--k-min
		Min k
		[Default: 1]
		
	--k-max
		Max k. Integer to specify a vaue or '' for infinity.
		[Default: '', i.e. infinity].	
		
	--max-GC
		Max GC content per motif (0.0 to 1.0)
		[Default: 1.0]
		
	--min-GC
		Min GC content per motif (0.0 to 1.0)	
		[Default: 1.0]
	
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
		[Default: 'bedtools merge']
	
	--bedtools-merge-opt = string
		A string with optional arguments that will 
		be passed to bedtools closest. 
		[Default: '']
		
	--print
		print out info while running.
	
	--help
		Print this man page.
		
# Examples:

	perl ./scripts/motif_search_vx-x.pl \
	-v test1 \
	-o ./MyData/Projects/GnNmSearch/2016m05d29-test \
	-f ./MyData/Projects/GnNmSearch/MotifSearch/ReferenceGenome/mm9/mm9.fa \
	--n-min 4 \
	--n-max '' \
	--m-min 1 \
	--m-max 4 \
	--k-min 4 \
	--k-max '' \
	--X-fwd '[G|g]' \
	--X-rev '[C|c]' \
	--print

