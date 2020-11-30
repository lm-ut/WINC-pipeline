### WINC-pipeline

### Introduction

WINC is a methodology that applies Local Ancestry concepts on ChromoPainter. WINC splits ChromoPainter copying vectors in windows with the same length and perform NNLS analyses on each windows. The Local Ancestry information will therefore be window-wise and not SNP-wise.  


All listed steps and scripts assume that the test dataset has been analyzed through ChromoPainter and the resulting copying vectors (expecially file.samples.out and file.recomrates) are available. Please run ChromoPainter setting your source populations as recipient (as you would do for the admixed population), in order to obtained copying vectors for them as well. You can set the source individuals as Recipient and as Donor in the same ChromoPainter run.   


All information about ChromoPainter can be found here: https://people.maths.bris.ac.uk/~madjl/finestructure-old/chromopainter_info.html  


### Splitting copying-vectors with ParseWindows script

The ParseWindows script allows to split ChromoPainter copying vectors into window-base copying vectors.  

ParseWindows is a python-based script and accepts the following parameters as arguments:  
  
  

*  -h, --help            show this help message and exit
*  --**label**           labelfile used for ChromoPainter, a three column file with ID information, population information and 0/1. 
*  --**samples**         ChromoPainter samples.out file 
*  --**genmap**          genetic map in bim format, a four columns file with Chr, SNP-ID, Position in cM, Base-pair coordinate
*  --**out**             output suffix  

Window length is currently embedded into the script and set as 500000 (500 kilo bases). To modify it please edit line 28 of the script and set the desired window length. This feature will be edited soon.  
  
  
    
Example of label file  
  
  

| Admixed_IND1 | Admixed | 1 |
|:------------:|:-------:|:-:|
| Admixed_IND2 | Admixed | 1 |
| Source1_IND1 | Source1 | 1 |
| Source2_IND1 | Source2 | 1 |
| INDn | Popn | 0 |  

                      
Example of command line  

`python ParseWindows.py --label ChromoPainter_DonorFile --samples ChromoPainter_IND1.samples.out --genmap BimFormat.GeneticMap --out GenWindows.out`  
  
    
    
#### Output

ParseWindows returns window-based copying vector per each file.sample.out analysed, where individuals are separated in two haplotypes (A and B).
The output files contain one row per window with the chunks donated to the individual in rows from populations in columns.  


Example 

| IND1 | StartingWindow | EndingWindow | DonorPop1 | DonorPop2 | DonorPopN |
|:----:|:--------------:|:------------:|:---------:|:---------:|:---------:|
| IND1 | 0 | 500000 | 0 | 0.40 | 0.85 | 2.7 |
| IND1 | 500000 | 1000000 | 0.15 | 0 | 1.35 |
