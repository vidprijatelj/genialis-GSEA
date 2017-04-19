# Gene Set Enrichment Analysis test code

Create and implement GSEA (Gene Set Enrichment Analysis)

### Prerequisites

Python 3.6.1; Pandas Library (0.19.0); Numpy Library (1.12.0);

Certain methods may not work outside the Linux/GNU ecosystem. Modifiy the file(s) accordingly

## Getting started

Copy the repository to your preferred location. 

All files should be in the same directory ("leukemia.txt" and "pathways.txt" included).

File "pathways.txt_clean.txt" can easily be ignored and/or deleted.

### How to run

Run file "gsea.py" (either "chmod 755 ./gsea.py" or by calling in "python3 gsea.py"

**PLEASE NOTE** - This program may take a while to complete. It takes ~1.5 seconds for each permutation. With 1000 permutations this may lead up to times of 25+ minutes.

Program ignores all the gene sets/pathways with "OBSOLETE" in their description or with the number of genes in the gene set of less than 15.

Wait for 1000 iterations to be completed.

Midway script outputs a file called "ES_table.csv" with all the calculated ES(S) and ES(S, pi) scores for all the included gene sets/pathways.

Ouput file is called "GSEA_statistics_$date-$time.csv". It is tab delimited and contains pathways and their NES(S) scores, p values and q values. Values are sorted by NES(S) in a descending manner (max -> min)

All NES(S) scores are bundled together (both positive and negative ones). 

Formatting was of my own choice and I did not follow the ones the authors used (Page 10: http://www.pnas.org/content/suppl/2005/10/14/0506580102.DC2/06580SuppText.pdf)
