import sys,gzip,pysam
from pysam import qualities_to_qualitystring as decode

workdir = sys.argv[1]

#infile = gzip.open(workdir+"/6_alignment.sam.gz","rb")
infile = pysam.AlignmentFile(workdir + "/6_alignment.bam","rb")
#outfile1 = gzip.open(workdir+"/7_alignment.sam.gz","wb")
outfile1 = pysam.AlignmentFile(workdir + "/7_alignment.bam","wb",template=infile)
outfile2 = gzip.open(workdir+"/8_rotated.fastq.gz","wb")

def Rotate(read):
	i = 0
	while i < read.query_length:
		outfile2.write("@" + read.query_name + "\n")
		outfile2.write(read.query_sequence[-i:] + read.query_sequence[:-i] + "\n")
		outfile2.write("+" + "\n")
		outfile2.write(decode(read.query_qualities[-i:]) + decode(read.query_qualities[:-i]) + "\n")
		i += 1

for read in infile.fetch(until_eof=True):
#for line in infile:
#	line1 = line
#	line = line.split()
	
	#make every possible rotation of reads with gaps and more than one clipped sequence
	cigar=read.cigarstring
	if read.is_unmapped:
		Rotate(read)
	elif cigar.count("D")> 0 or cigar.count("I") > 0 or cigar.count("S") > 1:
#	if line[5].count("D") > 0 or line[5].count("I") > 0 or line[5].count("S") > 1:
		Rotate(read)
	
	elif cigar.count("S") == 0:
#		xm = line[13].split(":") #xm = number of mismatches
		
		#save perfectly matched, ungapped alignments
		if read.get_tag("NM") == 0:
			outfile1.write(read)
		
		#make every possible rotation of reads with no gaps but have mismatches. Coincidentally matched bases near the
		#ends may suppress sequence clipping
		else:
			Rotate(read)
	
	#make every possible rotation of reads that still have clipped bases following rearrangement
	elif cigar.count("S") == 1:
		Rotate(read)

infile.close()
outfile1.close()
outfile2.close()
