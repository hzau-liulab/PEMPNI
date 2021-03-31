# PEMPNI
PEMPNI (Predictors for Effects of Mutations on PNIs) is a computational method that could effectively predict the effects of missense mutations on protein-DNA and protein-RNA interactions. 
This algorithm includes an energy module and a nonenergy module. The former is based on the novel geometric partition-based energy features, and the latter is based on the structural and sequence features. The integration of of multifaceted information could generate more robust predictions.

# Folder description
The ./bin directory contains programs to predict the two types of mutation effects.<br>
The ./datasets directory contains the mutational information of samples in various datasets.<br>
The ./WT_PDB directory contains the PDB files of wild-type structures used in this project.<br>
The ./MT_PDB directory contains the PDB files of mutant structures used in this project.<br>
The ./scripts directory contains auxiliary scripts used in this project.<br>
The ./features directory contains feature files used in this project.<br>
The ./examples directory contains two examples for protein-DNA/RNA complexes.<br>
The ./software directory contains software lists used in this project.<br>

# Installation
## Dependencies
Perl ( >= 5.0 ) and  pdl module<br>
Python ( >= 3.0 ) and scikit-learn, numpy, scipy, math and re modules<br>
GCC ( >= 5.3.0 )<br>

## Download
git clone https://github.com/hzau-liulab/PEMPNI.git

## Installation of third-party software 
To use PEMPNI, you need to conduct the installation of some third-party software, including Amber, Modeller, NACCESS, HBPLUS, ENDES and BLAST.<br>
1.To avoid copyright issues, we only provide the official website of third-party software. Please download and install them by yourself. It is recommended to install all the software in the same directory for modifying the paths.<br>
Amber: http://ambermd.org/index.php<br>
Modeller: https://salilab.org/modeller/<br>
NACCESS: http://wolf.bms.umist.ac.uk/naccess/<br>
HBPLUS: https://www.ebi.ac.uk/thornton-srv/software/HBPLUS/<br>
ENDES: https://sparks-lab.org/publication/ (Paper ID: 098)<br>
BLAST: ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/<br>
Nonredundant database: ftp://ftp.ncbi.nlm.nih.gov/blast/db/<br>

2.Related settings of Amber software<br>
Modify the first three lines of ./scripts/GB1/tleap.in and the second line of ./scripts/GB1/amber.sh according to the installation path of Amber software. Conduct the same settings for ./scripts/GB2/tleap.in and ./scripts/GB2/amber.sh.

3.Related settings of other software<br>
Other software should be installed in the same directory and obtain executable permissions. You can refer to the commands as follows.<br>
cd ./software<br>
chmod 777 ./hbplus/hbplus<br>
chmod 777 ./modeller9.25/bin/mod9.25<br>
chmod 777 ./naccess/naccess<br>
chmod 777 ./ncbi-blast-2.11.0+/bin/psiblast<br>
chmod 777 ./ncbi-blast-2.11.0+/bin/blastpgp<br>
chmod 777 ./enrich/processpdb<br>
chmod 777 ./pinup/pinup<br>
Note: Because processpdb and pinup use the blastpgp and nonredundant database, you need to manually modify the sixth line of these two files correspondingly.<br>

# How to use PEMPNI?
## Path settings
You will see protein-DNA.pl and protein-RNA.pl scripts in the ./bin/ directory, and you need to check and modify the paths related to specific software and databases involved in these scripts. For different types of complexes, you need to select the corresponding script.

## Input file
The input file is the original structure downloaded from the PDB database, and the scripts will automatically process the raw data. You should ensure that the structural file exists in the ./raw_PDB directory.

## Run PEMPNI
Please go to the ./bin/ directory, and then run the following command including the information of mutation sites. To run the script, you need to set three parameters, such as the PDB name, chain name and mutation information.
For example<br>
cd ./bin<br>
perl ./protein_DNA.pl 1AAY A D120A >>./log.txt<br>
or<br>
perl ./protein_RNA.pl 1FEU A D87E >>./log.txt<br>

# Results from PEMPNI
You will obtain all the resulting files in the  ./output/ directory, and the result.txt file contains the prediction results as follows.<br>
PDB ID:1AAY<br>
Chian ID:A<br>
Position:D120A<br>
Predicted affinity change:0.089<br>
Predicted score:0.122<br>
Significant decrease:NO<br>

# Interpretation of results
Predicted affinity changes are typically in the range of -3 to 5, and positive and negative values correspond to destabilizing and stabilizing effects, respectively. Predicted scores are in the range of 0 to 1, and if this score is greater than 0.47 for MPDs or 0.32 for MPRs, the submitted mutation could significantly decrease binding affinity.

# Help and Support
## Contact
If you have any questions or suggestions about PEMPNI, please contact us by email: jiangyao@webmail.hzau.edu.cn or liurong116@mail.hzau.edu.cn
The scripts of this project are free for academic use. For commercial use, please contact with the authors.
