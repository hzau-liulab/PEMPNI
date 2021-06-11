#! /user/bin/sh
path="/home/yaojiang/amber18";
$path/bin/tleap -f tleap.in >>log.txt;
$path/bin/sander -O -i min1.in -o min1.out -p com_solvated.prmtop -c com_solvated.inpcrd -r min1.rst -ref com_solvated.inpcrd;
$path/bin/sander -O -i min2.in -o min2.out -p com_solvated.prmtop -c min1.rst -r min2.rst -ref min1.rst;
$path/bin/sander -O -i min3.in -o min3.out -p com_solvated.prmtop -c min2.rst -r min3.rst -ref min2.rst;
$path/bin/cpptraj -i cpptraj.in >>log.txt;
$path/bin/MMPBSA.py -O -i mmpbsa_residue_pair.in -o Final_result_mmpbsa_GB2_pair.dat -do FINAL_DECOMP_MMPBSA_GB2_pair.dat -sp com_solvated.prmtop -cp com.prmtop -rp pro.prmtop -lp dna.prmtop -y last.crd > mmpbsa.log;
