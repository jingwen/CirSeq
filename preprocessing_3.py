import sys,gzip,pysam

workdir = sys.argv[1]

infiles = ["/9_alignment.bam","/10_alignment.bam"]
#outfile = gzip.open(workdir + "/11_alignment.sam.gz","wb")

infile = pysam.AlignmentFile(workdir + infiles[0],"rb")
outfile = pysam.AlignmentFile(workdir + "/11_alignment.bam","wb",template=infile)

#for file in infiles:
def CompareAS(infile,outfile):	
	#running lists
	Alignment = []
	AlignmentScore = []
	
	#current read
	SequenceID = ""
	
	for read in infile.fetch(until_eof=True):
#		RawLine = line
#		line = line.split()
		
		#add alignments of the same read to the running list
		if read.query_name == SequenceID:
			#add only alignments lacking gaps and clipped bases
			cigar=read.cigarstring
			if cigar.count("D") == 0 and cigar.count("I") == 0 and cigar.count("S") == 0:
				AS = read.get_tag("AS") 
				AlignmentScore.append(AS)
				Alignment.append(read)
				
		else:
			if len(Alignment) > 0:
				#write alignment with the best alignment score
				outfile.write(Alignment[AlignmentScore.index(max(AlignmentScore))])
			
			SequenceID = read.query_name
			
			Alignment = []
			AlignmentScore = []
			cigar=read.cigarstring
			if cigar.count("D") == 0 and cigar.count("I") == 0 and cigar.count("S") == 0:
				AS = read.get_tag("AS")
				AlignmentScore.append(AS)
				Alignment.append(read)



	if len(Alignment) > 0:
		#write alignment with the best alignment score
		outfile.write(Alignment[AlignmentScore.index(max(AlignmentScore))])

CompareAS(infile,outfile)
infile.close()
infile = pysam.AlignmentFile(workdir + infiles[1],"rb")
CompareAS(infile,outfile)
infile.close()
outfile.close()
