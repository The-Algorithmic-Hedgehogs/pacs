import numpy as np
import sys, math, os, csv
import argparse
import scipy.stats as stats

#this method creates the dictionaries that map guide sequence to number to reads and guide sequence to gene name
def createDictionaries(input_file):
	countDict = {}; #dictionary to keep track of the read counts for each guide
	guideGeneDict = {}; #dictionary to keep track of which gene each guide belongs to
	geneCountDict = {}; #dictionary to keep track of the read counts for a gene
	geneToGuideDict = {}; #dictionary to keep track of the sgRNAs for a gene
	# open library of guide sequences (.csv file) and create dictionary of read counts for each guide where count is initialized to 0
	try:
		with open(input_file) as csvfile:
			reader = csv.reader(csvfile, delimiter=',');
			next(reader); #skip first line
			for row in reader:
				guideSeq = row[7];
				countDict[guideSeq] = 0; #initialize sequence read count in countDict
				gene = row[2];
				guideGeneDict[guideSeq] = gene; #create guide->gene mapping in guideGeneDict
				geneCountDict[gene] = 0; #initialize read count for each gene
				#create gene->list of sgRNAs in geneToGuideDict
				if geneToGuideDict.get(gene) == None:
					geneToGuideDict[gene] = [guideSeq];
				else:
					(geneToGuideDict[gene]).append(guideSeq);
	except:
		print('could not open', input_file);
	return countDict, guideGeneDict, geneCountDict, geneToGuideDict;

#function to compute the edit distance between strings A and B using dynamic programming
# A Dynamic Programming based Python program for edit 
# distance problem 
def editDistDP(str1, str2, max):
	m = len(str1);
	n = len(str2);
	# Create a table to store results of subproblems
	dp = np.zeros( (m+1, n+1) );
	
	#fill in dp in top down manner
	for i in range(m+1):
		for j in range(n+1):
			# If first string is empty, only option is to insert all characters of second string
			if i == 0:
				dp[i][j] = j;
			# If second string is empty, only option is to remove all characters of second string
			elif j == 0:
				dp[i][j] = i; #minimum number of operations needed is i
			# If last characters are same, ignore last char and recur for remaining string
			elif str1[i-1] == str2[j-1]:
				dp[i][j] = dp[i-1][j-1];
			# If last character are different, consider all possibilities and find minimum
			else:
				dp[i][j] = 1 + min(dp[i][j-1], #insert
					dp[i-1][j], #remove
					dp[i-1][j-1]); #replace
		if dp[i][i] > max:
			return max + 1;

	return dp[m][n];
  

GUIDE_START = 0 #start index of guide sequence
GUIDE_END = 20 #end index of guide sequence

#matches the DNA sequencing reads to the guides. Returns the number of perfectly matched reads
def count_spacers(fastq_file, countDict, guideGeneDict, geneCountDict): 
	"""
	creates a dictionary with guide counts from fastq_file, writes to output_file
	fastq_file: forward read fastq file
	countDict: 
	dictionary: guide sequence as key, guide count as entry
	"""

	num_reads = 0 #total number of reads processed
	perfect_matches = 0 # guides with perfect match to library
	ed_matches = 0;

	# open fastq file
	try:
		ifile = open(fastq_file);
	except:
		print("could not find fastq file");
		return;

	# process reads in fastq file
	for seq in ifile: #contains the seq and Qscore etc.
		num_reads += 1;
		read_sequence = seq.strip(); #get rid of \n at the end
		for i in range(6): #go through first 6 bases to see if a guide can be located in the sequencing read
			start_idx = GUIDE_START + i; #get start index to look for guide
			end_idx = GUIDE_END + i; #get end index of sequence
			guide = read_sequence[start_idx : end_idx] #get possible sgRNA sequence from read
			if guide in countDict.keys():
				countDict[guide] += 1; #add one to count for guide
				perfect_matches += 1;
				gene = guideGeneDict[guide]; #find gene associated with guide
				geneCountDict[gene] += 1; #add one to count for gene
				break; #break from inner loop

			if i == 5: #no perfect match could be made so calculate edit distance
				guides = list(countDict.keys()); #get list of guides
				seq_to_compare = read_sequence[GUIDE_START : GUIDE_END];
				eds = np.zeros( len(guides) );
				for i, guide in enumerate(guides):
					eds[i] = editDistDP(seq_to_compare, guide, 3);

				guide_idx = np.argmin(eds); #get the index of the best matching guide
				matching_guide = guides[guide_idx]; #get best matching guide
				countDict[matching_guide] += 1; #add one to count for guide
				ed_matches += 1;
				gene = guideGeneDict[matching_guide]; #find gene associated with guide
				geneCountDict[gene] += 1;


	return perfect_matches, ed_matches;


#class to hold information about gene
class Gene:
	def __init__(self, name, unsortedReads, unsortedTotal, sortedReads, sortedTotal):
		self.name = name; #keeps track of gene name
		self.unsortedReads = unsortedReads;
		self.unsortedTotal = unsortedTotal;
		#sorted stats
		self.sortedReads = sortedReads;
		self.sortedTotal = sortedTotal;
		#p-values
		self.enrichPVal = 0;


#extract gene names and read counts from unsorted and sorted populations, and perform Fisher's exact test to determine enrichment p-value
def calcGeneEnrich(unsortedGeneCountDict, unsortedTotMatches, sortedGeneCountDict, sortedTotMatches):
	geneStatsDict = {}; #new dictionary to hold Gene objects
	keys = unsortedGeneCountDict.keys();
	for key in keys:
		#make new gene object with the correct counts
		gene = Gene(key, unsortedGeneCountDict[key], unsortedTotMatches, sortedGeneCountDict[key], sortedTotMatches);
		geneStatsDict[key] = gene;
		oddsratio, pValue = stats.fisher_exact([[gene.unsortedReads, gene.unsortedTotal], [gene.sortedReads, gene.sortedTotal]]);
		gene.enrichPVal = pValue;
	

	return geneStatsDict;

#prints the enrichment p-values for each gene to a .csv file. The file is ranked from lowest FDR-corrected p-value to highest
def printGeneEnrichStats(geneStatsDict, output_file):
	genesList = list( geneStatsDict.values() );
	genesList.sort(key = lambda gene: gene.enrichPVal);
	outputName = str(output_file+"_gene_enrichment_calculation.csv");
	with open(outputName, 'w') as csvfile:
		mywriter = csv.writer(csvfile, delimiter=',');
		mywriter.writerow(["Gene Name",  "p-value", "Reads in Unsorted Population", "Total Reads in Unsorted Population", "Reads in Sorted Population", "Total Reads in Sorted Population"]);
		for gene in genesList:
			mywriter.writerow([gene.name, gene.enrichPVal, gene.unsortedReads, gene.unsortedTotal, gene.sortedReads, gene.sortedTotal]);
		print("Printed Gene Enrichment stats");
	


#main method that ties processes together
def main(argv):
	parser = argparse.ArgumentParser(description='Analyze sequencing data for sgRNA library distribution');
	parser.add_argument('-s', '--sorted', type=str, dest='sorted_fastq', help="fastq file for the sorted population", default=None);
	parser.add_argument('-u', '--unsorted', type=str, dest='unsorted_fastq', help='fastq file for unsorted population', default=None);
	parser.add_argument('-o', '--output', type=str, dest='output_file',
						help='output file name', default='library');
	parser.add_argument('-g', '--guides', type=str, dest='guides_file',
						help='input file name', required=True);

	args = parser.parse_args();

	#necessary dictionaries to keep track of guide/gene info
	countDict, guideGeneDict, geneCountDict, geneToGuideDict = createDictionaries(args.guides_file);
	print("created required dictionaries");
	
	# count reads per guides and per gene for unsorted population if the fastq file is provided
	unsortedGuideCountDict = dict(countDict);
	unsortedGeneCountDict = dict(geneCountDict);
	print("Counting reads for unsorted population");
	unsortedTotMatches, unsorted_ed_matches =  count_spacers(args.unsorted_fastq, unsortedGuideCountDict, guideGeneDict, unsortedGeneCountDict);
	unsorted_total_reads = unsortedTotMatches + unsorted_ed_matches; #add perfect matches and matches by ED to get total number of reads
	print("Perfect matches in unsorted fastq: {0}".format(unsortedTotMatches));
	print("ED matches in unsorted fastq: {0}".format(unsorted_ed_matches));

	#count reads per guides and per gene for sorted population if the fastq file is provided
	sortedGuideCountDict = dict(countDict);
	sortedGeneCountDict = dict(geneCountDict);
	print("Counting reads for sorted population");
	sortedTotMatches, sorted_ed_matches  = count_spacers(args.sorted_fastq, sortedGuideCountDict, guideGeneDict, sortedGeneCountDict);
	sorted_total_reads = sortedTotMatches + sorted_ed_matches;  #add perfect matches and matches by ED to get total number of reads
	print("Perfect matches in sorted fastq: {0}".format(sortedTotMatches));
	print("ED matches in sorted fastq: {0}".format(sorted_ed_matches));

	#look for genes enriched for reads
	geneStatsDict = calcGeneEnrich(unsortedGeneCountDict, unsorted_total_reads, sortedGeneCountDict, sorted_total_reads);
	printGeneEnrichStats(geneStatsDict, args.output_file);

	###END OF MAIN


#run main method
main(sys.argv[1:]);