from pyspark import SparkConf, SparkContext
import sys, math, os, csv
import argparse
import scipy.stats as stats

conf = SparkConf().setAppName('pacs');
sc = SparkContext(conf = conf);


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

	print(args.sorted_fastq);
	print(args.unsorted_fastq);
	print(args.output_file);
	print(args.guides_file);

#run main method
main(sys.argv[1:]);



# pattern = sys.argv[1]; #get the pattern

# RDDvar = sc.textFile("ratings_full.csv"); #create RDD from the Gutenberg textfile

#  #use filter to find the lines with the specified pattern
# match_lines = RDDvar.filter(lambda line: pattern in line.lower());
# match_lines = match_lines.collect(); #changes RDD to list

# #write the matching lines to an output file
# with open("{0}_grep_output.txt".format(pattern), 'w') as ofile:
# 	for line in match_lines:
# 		ofile.write(line+'\n');
# 	ofile.close();