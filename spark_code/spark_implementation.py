from pyspark import SparkConf, SparkContext
import sys, math, os, csv
import argparse
import numpy as np
import scipy.stats as stats

### GLOBAL VARIABLES
conf = SparkConf().setMaster('local').setAppName('pacs');
sc = SparkContext(conf = conf);
GUIDE_START = 0 #start index of guide sequence
GUIDE_END = 20 #end index of guide sequence



#this method creates a dictionary that maps guide sequence to gene name and a dictionary of the genes in our database
def createDictionaries(input_file):
	guideGeneDict = {}; #dictionary to keep track of which gene each guide belongs to
	genesDict = {}; #dictionary to store gene name {'gene_name':0} for easy referencing
	# open library of guide sequences (.csv file) and create dictionary mapping guide sequence to gen
	try:
		with open(input_file) as csvfile:
			reader = csv.reader(csvfile, delimiter=',');
			next(reader); #skip first line
			for row in reader:
				guideSeq = row[7];
				gene = row[2];
				guideGeneDict[guideSeq] = gene; #create guide->gene mapping in guideGeneDict
				genesDict[gene] = 0; #create entry for gene in 'genesDict'

	except:
		print('could not open', input_file);
	return guideGeneDict, genesDict;


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


#function that matches the DNA sequencing read to one of the guides in the 'guideGeneDict'. 
#Returns a tuple (guide_sequence, 1)
def map_sequence(sequence, guideGeneDict): 
	
	for i in range(6): #go through first 6 bases to see if a guide can be located in the sequencing read
		start_idx = GUIDE_START + i; #get start index to look for guide
		end_idx = GUIDE_END + i; #get end index of sequence
		poss_guide = sequence[start_idx : end_idx] #get possible sgRNA sequence from read
		
		if poss_guide in guideGeneDict.keys(): #see if trimmed sequence is in dictionary of guides
			return (poss_guide, 1);

		if i == 5: #no perfect match could be made so calculate edit distance
			guides = list(guideGeneDict.keys()); #get list of guides
			seq_to_compare = sequence[GUIDE_START : GUIDE_END];
			eds = np.zeros( len(guides) );
			for i, guide in enumerate(guides):
				eds[i] = editDistDP(seq_to_compare, guide, 3);

			guide_idx = np.argmin(eds); #get the index of the best matching guide
			matching_guide = guides[guide_idx]; #get best matching guide
			
			return (matching_guide, 1);

#function that takes in a tuple 'guide_count', which is (guide_seq, count), and maps guide_seq to its corresponding gene in 'guideGeneDict'
# Returns a tuple (gene, count)
def map_guide_to_gene(guide_count, guideGeneDict):
	guide = guide_count[0]; #get guide sequence from tuple
	count = guide_count[1]; #get count for guide from tuple
	gene = guideGeneDict[guide]; #get corresponding gene for guide

	return (gene, count);


#class to hold the information regarding each gene
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



#main method that ties processes together
def main(argv):
	parser = argparse.ArgumentParser(description='Analyze sequencing data for sgRNA library distribution');
	parser.add_argument('-s', '--sorted', type=str, dest='sorted_fastq', help="fastq file for the sorted population", default=None);
	parser.add_argument('-u', '--unsorted', type=str, dest='unsorted_fastq', help='fastq file for unsorted population', default=None);
	parser.add_argument('-o', '--output', type=str, dest='output_file',
						help='output file name', default='output_from_pacs');
	parser.add_argument('-g', '--guides', type=str, dest='guides_file',
						help='input file name', required=True);

	#get the arguments from command-line
	args = parser.parse_args();

	#necessary dictionaries to keep track of guide/gene info
	guideGeneDict = createDictionaries(args.guides_file);

	### PERFORM MAPPING and EDIT DISTANCE for control sequences
	unsorted_rdd = sc.textFile(args.unsorted_fastq); #create RDD for unsorted reads
	unsorted_total_reads = unsorted_rdd.count(); #get number of total sequences in the file
	#map the sequences to guides, so 'unsorted_rdd' contains tuples (guide_seq, 1)
	unsorted_rdd = unsorted_rdd.map(lambda seq: map_sequence(seq, guideGeneDict));
	unsorted_guide_counts = unsorted_rdd.reduceByKey(lambda a, b: a+b); #get total count for each guide
	#map guide counts to genes
	unsorted_gene_counts = unsorted_guide_counts.map(lambda guide_count: map_guide_to_gene(guide_count, guideGeneDict));
	unsorted_gene_counts = unsorted_gene_counts.reduceByKey(lambda a, b: a+b); #get the total count for each gene
	unsorted_gene_counts = unsorted_gene_counts.collect(); #get a list of tuples (gene, count)

	### PERFORM MAPPING and EDIT DISTANCE for experimental sequences
	sorted_rdd = sc.textFile(args.sorted_fastq); #create RDD for sorted reads
	sorted_total_reads = unsorted_rdd.count(); #get number of total sequences in the file
	#map the sequences to guides, so 'sorted_rdd' contains tuples (guide_seq, 1)
	sorted_rdd = sorted_rdd.map(lambda seq: map_sequence(seq, guideGeneDict));
	sorted_guide_counts = sorted_rdd.reduceByKey(lambda a, b: a+b); #get total count for each guide
	#map guide counts to genes
	sorted_gene_counts = sorted_guide_counts.map(lambda guide_count: map_guide_to_gene(guide_count, guideGeneDict));
	sorted_gene_counts = sorted_gene_counts.reduceByKey(lambda a, b: a+b); #get the total count for each gene
	sorted_gene_counts = sorted_gene_counts.collect(); #get a list of tuples (gene, count)

	### PERFORM GENE ENRICHMENT ANALYSIS
	unsorted_gene_counts_dict = dict(unsorted_gene_counts); #create dictionary out of tuples
	sorted_gene_counts_dict = dict(sorted_gene_counts); #create dictionary out of tuples
	gene_names
	#calculation gene enrichment for 
	geneStatsDict = calcGeneEnrich(unsorted_gene_counts_dict, unsorted_total_reads, sorted_gene_counts_dict, sorted_total_reads);



	#write the counts of each URL to an output file
	# totCount = 0;
	# with open(args.output_file, 'w') as ofile:
	# 	for tup in unsorted_gene_counts:
	# 		ofile.write(tup[0] + '\t' + str(tup[1]) +'\n');
	# 		totCount += tup[1];
	# ofile.close();

	


#run main method
main(sys.argv[1:]);















