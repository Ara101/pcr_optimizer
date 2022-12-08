# pcr_optimizer
A Python package to aid researchers in optimizing pcr protocols

Polymerase chain reaction (PCR) is a common molecular biology technique used to amplify DNA. Each enzyme has different optimal reaction conditions such as extension times and primer annealing temperatures depending on the properties of the target gene, the primers used to amplify it, and buffers present in the reaction. PCR enzymes are expensive reagents, so it is also important to consider cost before choosing a protocol. Additionally, PCRs can be ineffective due to user design errors. Primers must anneal at the correct sites on a gene and have a carefully calculated annealing temperature to ensure proper binding conditions. Checking for any errors before conducting PCR is also crucial to avoid wasting time and money due to a failed reaction. Considering all these variables can make designing an effective PCR experiment difficult. 

This tool optimizes the PCR protocol decision making process. First the gene and primer sequences are checked for errors and evaluated to ensure the primers are binding in the correct location on the target sequence. Next, enzyme amount, cost per reaction, annealing temperature, annealing time, extension time, and total PCR reaction time for two common PCR enzymes (iProof High-Fidelity Polymerase and Taq Polymerase) are reported to the user. The user can also designate a factor (time or cost) to focus their optimization. 

## Setup 

```
pip install pcr-optimizer==0.13
pip install biopython
```

## Defining PCR object 

The gene, forward primer, and reverse primer sequences are entered as strings in 5’-3’ format. Template type must also be defined (plasmid, lambda, BAC DNA, or genomic), or will be assigned "plasmid" by default.  

Example: 

``` 
from pcr_optimizer.pcrprotocoloptimizer import pcr

gene = “atggagacagacacactcctgctatgggtactgctgctctgggttccaggttccactggtgacacaagtttgtacaaaaaagttggcaccaagtcgatcctagatggccttgcagataccaccttccgcaccatcaccactgacctcctgtacgtgggctcaaatgacattcagtacgaagacatcaaaggtgacatggcatccaaattagggtacttcccacagaaattccctttaacttcctttaggggaagtcccttccaagagaagatgactgcgggagacaacccccagctagtcccagcagaccaggtgaacattacagaattttacaacaagtctctctcgtccttcaaggagaatgaggagaacatccagtgtggggagaacttcatggacatagagtgtttcatggtcctgaaccccagccagcagctggccattgcagtcctgtccctcacgctgggcaccttcacggtcctggagaacctcctggtgctgtgcgtcatcctccactcccgcagcctccgctgcaggccttcctaccacttcatcggcagcctggcggtggcagacctcctggggagtgtcatttttgtctacagcttcattgacttccacgtgttccaccgcaaagatagccgcaacgtgtttctgttcaaactgggtggggtcacggcctccttcactgcctccgtgggcagcctgttcctcacagccatcgacaggtacatatccattcacaggcccctggcctataagaggattgtcaccaggcccaaggccgtggtggcgttttgcctgatgtggaccatagccattgtgatcgccgtgctgcctctcctgggctggaactgcgagaaactgcaatctgtttgctcagacattttcccacacattgatgaaacctacctgatgttctggatcggggtcaccagcgtactgcttctgttcatcgtgtatgcgtacatgtatattctctggaaggctcacagccacgccgtccgcatgattcagcgtaccgacgcgctggacctggaggagggaggaaacgtctatatcaaggccgacaagcagaagaacggcatcaaggcgaacttctgcatccgccacaacatcgaggacggcggcgtgcagctcgcctaccactaccagcagaacacccccatcggcgacggccccgtgctgctgcccgacaaccactacctgagcgtgcagtccaaactttcgaaagaccccaacgagaagcgcgatcacatggtcctgctggagttcgtgaccgccgccgggatcactttcggcatggacgagctgtacaagggcggtaccggagggagcatggtgagaaagggcgaggagctgttcaccggggtggtgcccatcctggtcgagctggacggcgacgtaaacggccacaagttcagcgtgggcggcgagggtgagggcgatgccaccgttggcaagctgaccctgaagttcatctgcaccaccggcaagctgcccgtgccctggcccaccctcgtgaccaccctgacctacggcgtgcagtgcttcagccgctaccccgaccacatgaagcagcacgacttcttcaagtccgccatgcccgaaggctacatccaggagcgcaccatcttcttcaaggacgacggcaactacaagacccgcgccgaggtgaagttcgagggcgacaccctggtgaaccgcatcgagctgaagggcatcgacttcaaggaggacggcaacatcctggggcacaagctggagtacaacaccggagcagcagcacgctggcgcgggcggcgcatggacattaggttagccaagaccctggtcctgatcctggtggtgttgatcatctgctggggccctctgcttgcaatcatggtgtatgatgtctttgggaagatgaacaagctcattaagacggtgtttgcattctgcaccatgctctgcctgctgaactccaccgtgaaccccatcatctatgctctgaggagtaaggacctgcgacacgctttccggagcatgtttccctcttgtgaaggcactgcgcagcctctggataacagcatgggggactcggactgcctgcacaaacacgcaaacaatgcagccagtgttcacagggccgcagaaagctgcatcaagagcacggtcaagattgccaaggtaaccatgtctgtgtccacagacacgtctgccgaggctctg"
forward_primer = "atggagacagacacactcctgctatgg"
reverse_primer = "cagagcctcggcagacgtgt"

my_pcr = pcr(gene, forward_primer, reverse_primer, template_type = "plasmid")
```
## Checking PCR object for errors: 

The check() function will check the gene and primer sequences for non-base charactrs (anything not A/T/G/C) and ensure primers are binding in the correct location. 

```
my_pcr.check()
```
Any errors will be reported and the user will manually fix them. 

## Optimzing PCR for cost or time: 

The recommend() function will optimize pcr for a user-defined factor (either "time" or "cost") and return a Dataframe of outputs. 

Optimizing for cost: 
```
my_pcr.recommend(factor = "cost")
```
```IProof Analyzer
+--------------------+-------+------------+-------+------------+--------+-------------+--------+-------------+
|  Reaction Volume   | 20 uL | 20 uL Cost | 50 uL | 50 uL Cost | 100 uL | 100 uL Cost | 200 uL | 200 uL Cost |
+--------------------+-------+------------+-------+------------+--------+-------------+--------+-------------+
| enzyme amount/cost |  0.2  |    0.31    |  0.5  |    0.77    |  1.0   |     1.53    |  2.0   |     3.06    |
+--------------------+-------+------------+-------+------------+--------+-------------+--------+-------------+
Taq Analyzer
+--------------------+-------+------------+-------+------------+--------+-------------+--------+-------------+
|  Reaction Volume   | 20 uL | 20 uL Cost | 50 uL | 50 uL Cost | 100 uL | 100 uL Cost | 200 uL | 200 uL Cost |
+--------------------+-------+------------+-------+------------+--------+-------------+--------+-------------+
| enzyme amount/cost |  0.1  |    0.18    |  0.25 |    0.45    |  0.5   |     0.89    |  1.0   |     1.78    |
+--------------------+-------+------------+-------+------------+--------+-------------+--------+-------------+
```

Optimizing for time: 
```
my_pcr.recommend(factor = "time")
```
```
IProof Analyzer
+-----------------------------+------------------------------+
|       Reaction Factor       |            Result            |
+-----------------------------+------------------------------+
|    Annealing Temperature:   |    56.72 degrees Celcius     |
|       Annealing Time:       |  30 seconds or 0.5 minutes   |
|       Extention Time:       | 1.5 seconds or 0.025 minutes |
| Total PCR reaction time is: | 34.71 minutes or 0.58 hours  |
+-----------------------------+------------------------------+
Taq Analyzer
+-----------------------------+------------------------------+
|       Reaction Factor       |            Result            |
+-----------------------------+------------------------------+
|    Annealing Temperature:   |    56.72 degrees Celcius     |
|       Annealing Time:       |  60 seconds or 1.0 minutes   |
|       Extention Time:       | 1.5 seconds or 0.025 minutes |
| Total PCR reaction time is: | 58.88 minutes or 0.98 hours  |
+-----------------------------+------------------------------+
````
If the user does not define a factor, then both tables will be returned for both enzymes. 
If the user defines factor as "return all" the output will be a tuple of all the optimization results
