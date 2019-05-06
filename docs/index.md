# Abstract

[CRISPR](https://en.wikipedia.org/wiki/CRISPR) screens are powerful experimental tools used to screen entire genomes in search of genes responsible for [phenotypes](https://en.wikipedia.org/wiki/Phenotype) of interest. With this approach, a single experiment can generate several gigabytes of data that, with sequential implementation, can take many hours to process, limiting [the](https://en.wikipedia.org/wiki/The) depth of sequencing and amount of analysis that can feasibly be performed. Here, we apply principles of parallel computing and algorithm design to expedite the data processing and analysis pipeline significantly. In doing so, we create a framework that provides three important features: i) Facilitation of considerably deeper sequencing experiments through parallelized expedition ii) integration of a sequence distance analysis for improved screening results iii) An associated cost-performance analysis for achieving the desired computational power within given financial constraints. 

# Introduction

## Background Information

The advent of technological advancements such as high-throughput sequencing and genome engineering, along with the increase in available computational power, has allowed biologists to adopt experimental approaches that create millions, sometimes even billions of data points per experiment. One example of such an advancement is the CRISPR genetic screen, in which researchers can introduce mutations in all twenty thousand (or so) genes in the human genome in parallel to identify new genes that are involved in a particular biological or pathological process (1).

The fundamental technology for this experimental procedure, as its name might suggest, is CRISPR genome engineering. This approach, developed from a bacterial viral defence mechanism, requires two components; a Cas9 protein, which acts like a set of molecular scissors that cuts DNA, and an RNA molecule that serves as a guide, ensuring that the Cas9 protein cuts at the specifically desired location. Using these two components, genes can be cut, copied and pasted with incredible efficiency and accuracy - and it is this efficiency that makes CRISPR genetic screens feasible.

The first step of the screen is to genetically engineer your cells of interest to express the Cas9 'scissor' protein, and to produce a large population of these cells. Next, you introduce a large library of different guide RNA molecules, which will stochastically insert themselves into your population of cells (with the help of some helpful little technique known as lentivirus transfection). The idea here is that every cell in your population has the potential for any of its genes to be removed (thanks to the presence of Cas9) - and through stochasticity, different cells will recieve guides and thus will cut out different genes. The result is a population of cells that have all the genes in their genomes mutated in parallel. 

Once this is achieved, you select the cells that display a phenotype of interest (how exactly this is done is very dependent on the particulars of the individual experiment), and sequence the DNA to identify which genes are mutated in these cells of interest.

In our particular experimental problem, we are looking at the [Shh](https://en.wikipedia.org/wiki/Sonic_hedgehog) pathway, a critical component of mammalian development. To identify cells in which this pathway is affected, the activity of this pathway is linked to the expression of a fluorescent protein, and cells with altered fluorescence are isolated with a technique known as [FACS](https://en.wikipedia.org/wiki/Flow_cytometry#Cell_sorting_by_flow_cytometry). In doing this, cells that have recieved mutations that affect this signalling pathway are identified and isolated.

## Problem Description

The output of the biological experiment for a CRISPR genetic screen consists of two files of DNA sequences:<br>
1) one file contains the DNA sequences from the control population of cells, (those that do not express the phenotype of interest) which we will call the *control file*<br>
2) the second file contains the DNA sequences from cells that were selected for some phenotype of interest, which we will call the *experimental file*.

The goal is to identify the mutations that are enriched in the experimental file compared to the control file.

Each file generally contains ten to twelve million DNA sequences that need to be processed. Each DNA sequence must be matched to a 'gold-standard' [reference database of ~80,000 sequences](https://github.com/rohuba/PACS/blob/master/data/Brie_CRISPR_library_with_control_guides.csv) to determine its origin. This matching process involves a string-to-string comparison. Generally, about 75% of the DNA sequences in a file can be perfectly matched to the gold standard. If the DNA sequence cannot be matched perfectly to this reference database, then the edit distance between the DNA sequence and each of the 80,000 gold-standard sequences must be calculated. String-to-string edit distance calculation is a very slow process, and it takes ~36s to calculate 80,000 edit distances for a single DNA sequence. Generally, around two million input sequences do not perfectly match the database of sequences, so calculating the edit distance between each of these sequences and the full reference database would take as much as 20,000 hours (36sec x 2M sequences / 3600 sec/hr). We would like to calculate the edit distances for the sequences that do not match any of the gold-standard sequences perfectly because this allows us to extract more information from a labor and time-intensive biological experiment.

Thus, this is a “Big Compute” problem because we need to compute these edit distances, and we are not sure which DNA sequences require edit distance calculation. Parallelization of this part of the application is required to make the application run in a reasonable amount of time (on the order of hours). Speed is also crucial because biological researchers depend on these results to execute their next experiments, which are time-sensitive because they are handling living cells and tissues.

MaGECK is an open-source pipeline for analyzing CRISPR screens (2) based on mean-variance modeling of the read counts for every gene. This pipeline provides read counts for genes given input control and experimental files, however it does not include a edit distance feature to determine the most likely gold-standard sequence for a sequencing read that did not perfectly match a gold-standard sequence. Additionally, MaGECK is considered a black-box program with statistics that are not well understood by general biologists. Our application would provide a missing functionality that would allow researchers to extract more information from their screens and give researchers a tool that can be clearly explained.


## Existing Pipeline Analysis
The existing pipeline takes in a *control file* and *experiment file* generated from the DNA sequencing step of the genetic screen. It also takes in a file containing the "database" of 80,000 gold-standard sequences. First, every DNA sequence in the *control file* is mapped to one of the gold-standard sequences; if the sequence matches one of the gold-standard sequences, then the match count for gold-standard sequence is incremented and the match count for corresponding gene is incremented. If no perfect match can be determined, then the edit distance between the control DNA sequence and every gold-standard sequence is calculated. The gold-standard DNA sequence with the lowest edit distance is found, and its match count and the match count for the coressponding gene are incremented. Once all the DNA sequences in the *control file* have been matched, this process is repeated for the *experimental file*. The match counts for genes in the *control file* and the *experiment file* are aggregated, and a Fisher's exact test is run on the match counts from both files to determine if there is a significant change in the number of matched sequences for each gene in the *experimental file* when compared to the *control file*. The output of the pipeline consists of, for each gene, its Fisher's exact test p-value, the number of matching sequences in the *control file*, the total number of sequences in the *control file*, the number of matching sequences in the *experimental file*, and the total number of sequences in the *experimental file*. These stats are written to a CSV file.

The code for the exisiting sequential pipeline for analyzing the results of CRISPR genetic screens can be found [here](https://github.com/rohuba/PACS/blob/master/sequential_pipeline/sequential_analysis.py). The Python script `sequential_analysis.py` takes in multiple arguments through the use of command-line flags:<br>
    1. `-u` indicates the file following file path if for the *control file*<br>
    2. `-s` indicates the file following file path if for the *experimental file*<br>
    3. `-g` indicates the file following file path if for a CSV file containing the 80,000 gold-standard sequences<br>
    4. `-o` indicates the following string will be used as a prefix for the output file<br>
    
Command-line Example:<br>
`python sequential_analysis.py -u ../data/control_file_100_seqs.txt -s ../data/experimental_file_100_seqs.txt -g ../data/Brie_CRISPR_library_with_controls_guides.csv -o test_output`

The output of the script is a CSV file with six columns with the following information:
    1. Name of the gene<br>
    2. $p$-value from Fisher's exact test for enrichment<br>
    3. Number of DNA sequences from *control file* mapping to gold-standard sequences contained in this gene<br>
    4. Total number of DNA sequences in *control file*
    5. Number of DNA sequences from *experimental file* mapping to gold-standard sequences contained in this gene<br>
    6. Total number of DNA sequences in *experimental file*
    
# Project Design

## Sequential Code Profiling

The sequential code `sequential_analysis.py` was profiled using the `cProfile` Python package. The file was run with a control file of 100 sequences (*control_file_100_seqs.txt*) and an experimental file of 100 sequences (*experimental_file_100_seqs.txt*). Each of these input files contained 75 sequencing reads that could be perfectly matched to the database of 80,000 guide sequences and 25 sequencing reads that needed an edit distance calculation. This breakdown was representative of the proportion of sequencing reads in the full input files; ~25% of sequencing reads cannot be perfectly matched to one of the 80,000 guide sequences. The exact command that was run was: `python -m cProfile -o 100_seq_stats.profile sequential_analysis.py -g ../data/Brie_CRISPR_library_with_controls_FOR_ANALYSIS.csv -u ../data/Genome-Pos-3T3-Unsorted_100_seqs.txt -s ../data/Genome-Pos-3T3-Bot10_100_seqs.txt -o cProfile_test_output`. This code was run on a Macbook Pro, with a 2.2 GHz Intel Core i7 processor with 6 cores. The profiling information was saved in a file called *100_seq_stats.profile*. The following results are from the `pstats` package.
```python
import pstats

p = pstats.Stats('100_seq_stats.profile'); #read in profiling stats
p.strip_dirs(); #remove the extraneous path from all the module names

#sort according to time spent within each function, and then print the statistics for the top 20 functions. 
p.sort_stats('time').print_stats(20)

```
	Mon Apr 29 16:08:45 2019    100_seq_stats.profile
	         
	         350132307 function calls (350126604 primitive calls) in 537.388 seconds

	   Ordered by: internal time
	   List reduced from 1999 to 20 due to restriction <20>

	   ncalls  tottime  percall  cumtime  percall filename:lineno(function)
	  3981650  462.315    0.000  534.561    0.000 count_spacers_with_ED.py:35(editDistDP)
	333437825   69.985    0.000   69.985    0.000 {built-in method builtins.min}
	  3982102    1.819    0.000    1.819    0.000 {built-in method numpy.zeros}
	        2    1.414    0.707  535.981  267.990 count_spacers_with_ED.py:69(count_spacers)
	7992171/7992125    0.447    0.000    0.447    0.000 {built-in method builtins.len}
	        1    0.215    0.215    0.232    0.232 count_spacers_with_ED.py:7(createDictionaries)
	    81/79    0.153    0.002    0.156    0.002 {built-in method _imp.create_dynamic}
	    20674    0.115    0.000    0.365    0.000 stats.py:3055(fisher_exact)
	      348    0.100    0.000    0.100    0.000 {method 'read' of '_io.FileIO' objects}
	    41950    0.092    0.000    0.092    0.000 {method 'reduce' of 'numpy.ufunc' objects}
	        1    0.075    0.075    0.075    0.075 {method 'dot' of 'numpy.ndarray' objects}
	      348    0.056    0.000    0.155    0.000 <frozen importlib._bootstrap_external>:830(get_data)
	    28062    0.052    0.000    0.052    0.000 {built-in method numpy.array}
	     1604    0.035    0.000    0.035    0.000 {built-in method posix.stat}
	      348    0.034    0.000    0.034    0.000 {built-in method marshal.loads}
	    593/1    0.033    0.000  537.389  537.389 {built-in method builtins.exec}
	        1    0.033    0.033    0.405    0.405 count_spacers_with_ED.py:135(calcGeneEnrich)
	    81/65    0.024    0.000    0.077    0.001 {built-in method _imp.exec_dynamic}
	    21078    0.021    0.000    0.071    0.000 fromnumeric.py:69(_wrapreduction)
	    20675    0.020    0.000    0.020    0.000 {method 'writerow' of '_csv.writer' objects}

The majority of runtime is spent with `editDistDP` function. 534 of the 537 seconds, which accounts for 99.4% of the runtime, are spent calculating the edit distance between 50 sequencing reads and 80,000 guides. Generally, the input files contain ~10M sequencing reads, and about 25% of the sequences cannot be matched perfectly to one of the 80,000 guides. Thus for two input files of ~10M sequencing reads (~20M reads total), there are ~4-5M sequencing reads for which the edit distance calculations must be performed. If this code was run sequentially, this would require 10,000 hours of runtime. Therefore, we need to parallelize this portion of the code.

The edit distance calculation is currently nested within the function `count_spacers`, which matches each sequencing read from the input files to one of the 80,000 guides. For 200 sequencing reads provided as input, 1.4 seconds are spent performing the matching. This is only 0.007 seconds per sequencing read (using the 1.4 seconds from the *tottime* column since the *cumtime* takes into account the edit distance calculation). This number grows large if we have 20M sequencing reads - it would take $\dfrac{0.007\text{seconds/read} \cdot 20\text{M reads}}{3600\text{seconds/hour}} = 39\text{hours}$. Thus, the entire matching process of our workflow needs to be parallelized.

We want to parallelize this matching process by using a Spark cluster to have access to as many cores as possible to perform both the matching process and edit distance calculation (if needed). We will partition each input file into many tasks, and each task will run on a single core of the Spark cluster. A single core will perform both the matching process and edit distance for the sequencing reads in a partition. From what we have determined, there is not an easy way to parallelize the edit distance calculation algorithm itself. However, for a given sequencing read, we should be able to parallelize the 80,000 edit distance calculations that need to be performed between the sequencing read and the guides by using Python multi-threading.

## Overheads 

Since we do not know for which sequences we will need to perform edit distance calculations, load-balancing is the main overhead we anticipate dealing with because we do not want one or two cores slowed down with having to compute too many edit distance calculations. We would like to spread the number of edit distance calculations out evenly between cores by tuning the number of Spark tasks. It may be good to shuffle the order of the sequencing reads because sometimes many sequences that require edit distance calculations are adjacent to each other in the input file.

If we try to use a GPU to perform the 80,000 edit distance calculations in parallel, memory-transfer (input/output) to the GPU would also be an overhead. For a single sequencing read, multiple transfers would need to be performed as we would not be able to perform the 80,000 calculations in parallel, since we are limited by the number of cores on the GPU. Currently, we do not have a good way of mitigating this GPU overhead.

## Scaling

The sequence matching and edit distance portion of our code accounts for 99.4% of the runtime in our small example. With larger problem sizes, this percentage should only increase because the number of operations performed after the sequence matching and edit distance section is constant.

Amdahl's Law states that potential program speedup $S_t$ is defined by the fraction of code $c$ that can be parallelized, according to the formula
$$
S_t = \dfrac{1}{(1-c)+\frac{c}{p}}
$$
where $p$ is the number of processors/cores. In our case, $c=0.994$ and the table below shows the speed-ups for $2$, $4$, $8$, $64$, and $128$ processors/cores:

|processors|speed-up|
|----------|--------|
|2|1.98x|
|4|3.93x|
|8|7.68x|
|64|46.44x|
|128|72.64x|

Thus, our strong-scaling is almost linear when $p$ is small, but we observe that this begins to break down because we only get 73x speed-up if we were to use 128 processors/cores.

Gustafson's Law states larger systems should be used to solve larger problems because there should ideally be a fixed amount of parallel work per processor. The speed-up $S_t$ is calculated by
$$
S_t = 1 - c +c\cdot p
$$
where $p$ is the number of processors/cores. In our case, $c=0.994$ and the table below shows the speed-ups for $2$, $4$, $8$, $64$, and $128$ processors/cores:

|processors|speed-up|
|----------|--------|
|2|1.994x|
|4|3.98x|
|8|7.96x|
|64|63.62x|
|128|127.23x|

Thus, we almost achieve perfect weak-scaling because we can split up larger problem-sizes (which would be larger input files in our case) over more processors to achieve about the same runtime.


## Data

TO DO

## Pipeline

TO DO

## Speedup Algorithm

TO DO


## Infrastructure

TO DO

* * *

# Usage 

TO DO

* * *

# Results

## Performance Evaluation

TO DO

## Optimizations and Overheads

TO DO

* * *

# Discussion

TO DO

# Future Work

3. Analyze control and experimental files on separate clusters and possibly use some message passing

# References

1. Pusapati GV, Kong JH, Patel BB, Krishnan A, Sagner A,
Kinnebrew M, Briscoe J, Aravind L, Rohatgi R: CRISPR screens
uncover genes that regulate target cell sensitivity to the
morphogen sonic hedgehog. Dev Cell 2018, 44:113-129 e118.<br>
2. Li, W., Xu, H., Xiao, T., Cong, L., Love, M.I., Zhang, F., Irizarry, R.A., Liu, J.S., Brown, M., and Liu, X.S. (2014). MAGeCK enables robust identification of essential genes from genome-scale CRISPR/Cas9 knockout screens. Genome Biol. 15, 554.<br>