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
*  --**out**             output name  

Window length is currently embedded into the script and set as 500000 (500 kilo bases). To modify it please edit line 28 of the script and set the desired window length. This feature will be edited soon.  
  
  
    
Example of label file  
  
  

| Admixed_IND1 | Admixed | 1 |
|:------------:|:-------:|:-:|
| Admixed_IND2 | Admixed | 1 |
| Source1_IND1 | Source1 | 1 |
| Source2_IND1 | Source2 | 1 |
| INDn | Popn | 0 |  

                      
Example of command line  

`python ParseWindows.py --label ChromoPainter_DonorFile --samples ChromoPainter_IND1.samples.out --genmap BimFormat.GeneticMap --out IND1GenWindows.out`  
  
    
    
#### Output

ParseWindows returns window-based copying vector per each individual analised in the file.sample.out. In each output file individuals are separated in two haplotypes (A and B).
The output files contain one row per window with the chunks donated to the individual in rows from populations in columns.  


Example 

| IND1 | StartingWindow | EndingWindow | DonorPop1 | DonorPop2 | DonorPopN |
|:----:|:--------------:|:------------:|:---------:|:---------:|:---------:|
| IND1 | 0 | 500000 | 0 | 0.40 | 0.85 | 2.7 |
| IND1 | 500000 | 1000000 | 0.15 | 0 | 1.35 |  

### Performing NNLS analyses with the addition of the Correlation-AssignmentScore matrix

NNLS_step.R is a R-based script that allows to perform NNLS on a window-based copying vector. The script applies the C-AS matrix as well (CASmatrix.txt), a reference grid to inform a priori on the accuracy to be expected by WINC for a given set of Ancestry Assignments and local diversification between sources.  

The script accepts a .csv file as argument, with the following parameters:  

* **S1**	             name of source 1 population
* **S2**	             name of source 2 population
* **S3**               OPTIONAL parameter with the name of source 3 population
* **admixedName**	     name of the admixed population
* **pathPainting**	   path to the admixed population copying vectors splitted by the ParseWindow 
* **donorDir**	       path to the source populations copying vectors splitted by the ParseWindow 
* **listSamplesFile**	 path to a 2 columns file,  with individuals' ID in the first colum and individuals' population infomation in the second one
* **suffixWindowsPainted** suffix used in the ParseWindows.py	--out argument

Example  

| S1 | S2 | admixedName | pathPainting | donorDir | listSamplesFile |
|:--:|:--:|:-----------:|:------------:|:--------:|:---------------:|
| CEU  |  ESN  |   ASW  |  ~/path1/  | ~/path2/ | ~/path2/SampleList.txt |


#### Output 

The NNLS_step script outputs a window-based Local Ancestry call per each individual (INDn.winc.tmp). The output files contain one row per window with the source assignments. The first two columns contain windows coordinates (Start and End). The following columns (S1 and S2 in the example, and S3 in case of a three-ways admixture) yield the proportion of the sources.  
Following the source columns information, we added the ancestry assignment calls obtained by applying either the C-AS matrix (columns named accuracy.N) or a simple threshold (colums named thresN). The C-AS matrix accuracy values used are: 0.8, 0.9, 0.95 and 0.99. The threshold values used are: 0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.91,0.92,0.93,0.94,0.95,0.96,0.97,0.98,0.99.  
The rho column contains the rho values obtained between sources' windows. The last column "Usable" checks if there is enough information in the window copying vector analysed.

The current output will be edited soon to be easier to handle.  

Example of output file  


| Start | End | S1 | S2 | NoThreshold | accuracy.8 | accuracy.N | thres0.55 | thresN | rho | Usable | 
|:-----:|:-----:|:-----:|:-----:|:-----:|:-----:|:-----:|:-----:|:-----:|:-----:|:-----:|
|   0   | 500000|1|0|0|0|0|0|0|-0.11|Y|
|500000|1000000|0.74|0.26|0|0|0|NA|NA|-0.12|Y|  

 
 
#### How to cite this work & get in touch

Molinaro et al. 2020 submitted for publication  

ludovica.molinaro@ut.ee
