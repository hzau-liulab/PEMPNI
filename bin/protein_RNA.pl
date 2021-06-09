use warnings;
use strict;
###@author jiang yao

#Examples:perl protein_RNA.pl 1FEU A D87E
#Notice:You need to set the correct working path and software installation path
##############################################################Set working directory#################################################
my $workdir="/home/yaojiang/protein-RNA/";             #Set to your own working directory 
my $software_dir="/var/www/html/PEMPNI/software/";     #Set the software installation directory
#
system "mkdir -p $workdir/output/data_process";
system "mkdir -p $workdir/output/energy_feature/GB2";
system "mkdir -p $workdir/output/distance";
system "mkdir -p $workdir/output/bio_feature/ASA";
system "mkdir -p $workdir/output/bio_feature/HB";
system "mkdir -p $workdir/output/bio_feature/JSD/pssm";
system "mkdir -p $workdir/output/bio_feature/JSD/cscores";
system "mkdir -p $workdir/output//bio_feature/ENDES";
#####################################################################################################################################
my %aa=('ALA'=>'A','ARG'=>'R','ASN'=>'N','ASP'=>'D','CYS'=>'C','GLN'=>'Q','GLU'=>'E','GLY'=>'G','HIS'=>'H','ILE'=>'I',
        'LEU'=>'L','LYS'=>'K','MET'=>'M','PHE'=>'F','PRO'=>'P','SER'=>'S','THR'=>'T','TRP'=>'W','TYR'=>'Y','VAL'=>'V');
my %reverse_aa=reverse %aa;

my $pdb_name = $ARGV[0];               #1FEU
my $chain_name = $ARGV[1];             #A
my $mutate_res = $ARGV[2];             #D87E
my $mutated_res=substr($mutate_res,-1,1);
my $mutate_restyp = $reverse_aa{$mutated_res}; #GLU
my $pdb_chain = $pdb_name.'_'.$chain_name;     #1FEU_A
my $mut_name="$pdb_chain.mut.$mutate_res";     #1FEU_A.mut.D87E
my $mutate_pos =$mutate_res;
   $mutate_pos =~s/[A-Z]//g;             #87
my $wild_res=substr($mutate_res,0,1);    #D
my $wlid_restyp=$reverse_aa{$wild_res};  #ASP

my @name;
push(@name,$pdb_chain);
push(@name,$mut_name);

system "cp $workdir/raw_PDB/$pdb_name.pdb $workdir/output/data_process";

my %AGCT=('A'=>1,'G'=>1,'C'=>1,'U'=>1);
########################################################Original structure processing######################################
open(IN,"$workdir/output/data_process/$pdb_name.pdb") or die "can not open this files\n"; 
open(OUT,">$workdir/output/data_process/WT_pro.pdb") or die "can not open $!\n";
open(OUTTTTT,">>$workdir/output/data_process/WT_dna.pdb") or die "can not open $!\n";
open(OUTT,">$workdir/output/data_process/pre_dna.pdb") or die "can not open $!\n";
open(OUTTT,">>$workdir/output/data_process/WT_com.pdb") or die "can not open $!\n";
open(OUTTTT,">>$workdir/output/data_process/$pdb_chain.pdb") or die "can not open $!\n";
while (my $line=<IN>) {
    if($line=~/^ATOM/){
    my $rna = substr($line,18,3);
    $rna =~s/\s+//g;
    my $atom = substr($line,77,1);
    $atom =~s/\s+//g;
    if ( exists($AGCT{$rna})) {print OUTT "$line";}
    elsif( $ARGV[1] eq substr($line,21,1) and $atom ne 'H') {print OUT "$line";}
     }
    }
close IN;
close OUT;
close OUTT;

my %hhh;
open(IN, "$workdir/output/data_process/WT_pro.pdb") or die "can not open $!\n";
my @propro=<IN>;
close IN;

open(IN, "$workdir/output/data_process/pre_dna.pdb") or die "can not open $!\n";
my @nulnul=<IN>;
close IN;

open(OOUT,">$workdir/output/data_process/pree_dna.pdb") or die "can not open $!\n";

foreach my $key1 (@nulnul){
    my $x = substr($key1,30,8);
    my $y = substr($key1,38,8);
    my $z = substr($key1,46,8);
    my $chain1 = substr($key1,21,1);
        foreach my $key2 (@propro){
           my $x1 = substr($key2,30,8);
           my $y1 = substr($key2,38,8);
           my $z1 = substr($key2,46,8);
           my $dis = sqrt(($x - $x1)**2+($y - $y1)**2+($z-$z1)**2);
                        if ($dis <= 100) {$hhh{$chain1}+=1;}
        }
}

my @kkk=keys %hhh;
foreach my $line (@nulnul){
    my $chain1 = substr($line,21,1);
    if (grep {$_ eq $chain1 } @kkk ) {
        print OOUT "$line";
    }   
}
close OOUT;

open(INN, "$workdir/output/data_process/pree_dna.pdb") or die "can not open $!\n";
my @preedna=<INN>;
close INN;

foreach my $chain (@kkk){
    my @dna_chain;
    foreach my $line (@preedna){
    my $chain1 = substr($line,21,1);
    if ($chain1 eq $chain) {
        push (@dna_chain,$line);
    }
    }
    my $num;
    foreach my $li (0..$#dna_chain){
        my $pp=substr($dna_chain[$li],13,3);
        $pp =~s/\s+//g;
         if ($pp eq "O5'") {
        $num=$li;
        last;    
    }
         }
    foreach my $i ($num..$#dna_chain){
        print OUTTTTT "$dna_chain[$i]";
    }
}
close OUTTTTT;

open(IN, "$workdir/output/data_process/WT_dna.pdb") or die "can not open $!\n";
my @nulll=<IN>;
close IN;

foreach my $line (@nulll){
    print OUTTT "$line";
    print OUTTTT "$line";
}
foreach my $line (@propro){
    print OUTTT "$line";
    print OUTTTT "$line";
}

close OUTTT;
close OUTTTT;
###########################################Check whether the file meets the amber calculation requirements#####################
my @flag;
for my $i (@kkk){
    my @nucleic;
    for my $j (@nulll){
        my $chai = substr($j,21,1);
        $chai =~s/\s+//g;
        if ($chai eq $i) {
           push (@nucleic,$j);
        }
    }
    my $pp=substr($nucleic[0],13,3);
    $pp =~s/\s+//g;
    if ($pp eq "O5'") {
     push (@flag,0);
}else{
     push (@flag,1);
}
}
my $sum=0;
for my $i (@flag){
    $sum+=$i;
}
if ($sum==0) {
    print "Your file meets Amber calculation requirements\n"
}else{
    print "Your file does not meet Amber calculation requirements, please reprocess\n";
}

##############################################################modeller constructs the mutation structure and processes##################################################
my $out = "$workdir/output/data_process/";
system "$software_dir/modeller9.25/bin/mod9.25 $workdir/scripts/mutate_model.py $out $mutate_res $pdb_chain $mutate_pos $mutate_restyp $chain_name >>./log.txt";

open(IN,"$workdir/output/data_process/$mut_name.pdb") or die "can not open this files\n"; 
open(OUT,">$workdir/output/data_process/MT_pro.pdb") or die "can not open $!\n";
open(OUTT,">$workdir/output/data_process/MT_dna.pdb") or die "can not open $!\n";
open(OUTTT,">$workdir/output/data_process/MT_com.pdb") or die "can not open $!\n";
open(OUTTTT,">$workdir/output/data_process/MT_$pdb_chain.pdb") or die "can not open $!\n";
while (my $line=<IN>) {
    if($line=~/^ATOM/){  
    my $rna = substr($line,18,3);
    $rna =~s/\s+//g;
    if ( exists($AGCT{$rna})) {print OUTT "$line"; print OUTTT "$line";print OUTTTT "$line";}
    elsif( $ARGV[1] eq substr($line,21,1)) {print OUT "$line";print OUTTT "$line";
                                           print OUTTTT "$line";}
     }
    }
close IN;
close OUT;
close OUTT;
close OUTTT;
close OUTTTT;

open(OUTTTTT,">>$workdir/output/data_process/DNA.pdb") or die "can not open $!\n";
my @nul;
open(INN,"$workdir/output/data_process/$mut_name.pdb") or die "can not open this files\n"; 
chomp(my @pdb = <INN>);
close INN;
my $numm;
for my $i (0..$#pdb){
    my $ress = substr($pdb[$i],17,3);
    if(exists $aa{$ress} && $pdb[$i]=~/^ATOM/){
      $numm = $i;
      last;
}
}

for my $i (6..$numm-2){
    push(@nul,$pdb[$i])
}

my %ter_num;
for my $i(0..$#nul){
    if ($nul[$i]=~/^TER/) {
       $ter_num{$i}=1;
       $ter_num{$i+1}=1;
       $ter_num{$i+2}=1;
       $ter_num{$i+3}=1;       
    }
}
for my $i(0..$#nul){
    if (exists $ter_num{$i}) {
       next;
    }else{
      print OUTTTTT "$nul[$i]\n";
    }
    
}
close OUTTTTT;

open(OU,">>$workdir/output/data_process/COM.pdb") or die "can not open $!\n";
open(IN,"$workdir/output/data_process/DNA.pdb") or die "can not open $!\n";
while (my $line = <IN>) {
    chomp $line;
    print OU "$line\n";
}
close IN;

open(IN,"$workdir/output/data_process/MT_pro.pdb") or die "can not open $!\n";
while (my $line = <IN>) {
    chomp $line;
    print OU "$line\n";
}
close IN;
close OU;

##############################################################Construct the fasta sequence of the protein##################################################
my @na;
push (@na,"MT_pro");
push (@na,"WT_pro");
for my $na (@na){
open(IN,"$workdir/output/data_process/$na.pdb") or die "can not open this files\n";
my %hash;
while (my $line=<IN>) {
   my $b=substr ($line,17,3);
   my $c=substr ($line,22,5);
   $c=~s/\s+//;
   $hash{$c}=$b;
}
close IN;

my @fasta;
for my $key (sort{$a<=>$b} keys %hash)
{
    my $value=&text($hash{$key});
    push(@fasta,$value);
}
open(OUT,">$workdir/output/data_process/$na.fas") or die "can not open this files\n";
print OUT ">"."$na\n";
foreach my $i (@fasta){
    print OUT "$i"
}
close OUT;

my $nn=1;
open(OUTT,">$workdir/output/data_process/$na"."_seq.txt") or die "can not open this files\n";
for my $key (sort{$a<=>$b} keys %hash)
{
    my $value=&text($hash{$key});
    print OUTT "$nn\t$value\t$key\n";
    $nn+=1;
}
close OUTT;

sub text {
my($input) = $_[0];
my %hash = (
'ALA' => 'A',
'VAL' => 'V',
'LEU' => 'L',
'ILE' => 'I',
'PRO' => 'P',
'TRP' => 'W',
'PHE' => 'F',
'MET' => 'M',
'GLY' => 'G',
'SER' => 'S',
'THR' => 'T',
'TYR' => 'Y',
'CYS' => 'C',
'ASN' => 'N',
'GLN' => 'Q',
'LYS' => 'K',
'ARG' => 'R',
'HIS' => 'H',
'ASP' => 'D',
'GLU' => 'E',
);
my $out=$hash{$input};
return $out;
}
}

#################################Calculate the file whose protein-nucleic acid atom distance is less than 5#################################################
for my $i (@name){
my @he_suan;
my @protein;

my @acid = keys %aa;
open(IN,"<","$workdir/output/data_process/$i.pdb") or die "$!";
while (my $line = <IN>) {
    if ($line =~/^ATOM/) {
    my $aa=substr($line,17,3);
    $aa =~s/\s+//g; 
        if (grep {$_ eq $aa}@acid) {
            push (@protein,$line)
           }else{
           push (@he_suan,$line)}
    }else{next;}     
}
close IN;

open(OUT,">","$workdir/output/distance/$i"."_pro_nul_dis_5.txt") or die "$!";
foreach my $key1 (@protein){
    my $key = substr($key1,23,3);
    $key=~s/\s+//;
    my $x = substr($key1,30,8);
    my $y = substr($key1,38,8);
    my $z = substr($key1,46,8);
    my $atom1 = substr($key1,13,3);
    my $res1 = $aa{substr($key1,17,3)};
    my $chain1 = substr($key1,21,1);
        foreach my $key2 (@he_suan){
           my $x1 = substr($key2,30,8);
           my $y1 = substr($key2,38,8);
           my $z1 = substr($key2,46,8);
           my $atom2 = substr($key2,13,3);
           my $base = substr($key2,17,3);
           my $numres2 = substr($key2,22,5);
           my $chain2 = substr($key2,21,1);
               $numres2 =~ s/\s+//g;
           my $dis = sqrt(($x - $x1)**2+($y - $y1)**2+($z-$z1)**2);
                        if ($dis <= 5) {
                            print OUT "$key\t$chain1\t$dis\t$res1\t$atom1\t$numres2\t$chain2\t$base\t$atom2\n";
                        }else{next;}
        }
    }
close OUT;
#################################Calculate files with protein-protein atom distance less than 4.5#################################################
open(OUT,">","$workdir/output/distance/$i"."_pro_pro_dis_4.5.txt") or die "$!";
foreach my $key1 (@protein){
    my $key = substr($key1,23,3);
    $key=~s/\s+//;
    my $x = substr($key1,30,8);
    my $y = substr($key1,38,8);
    my $z = substr($key1,46,8);
    my $atom1 = substr($key1,13,3);
    my $res1 = $aa{substr($key1,17,3)};
    my $chain1 = substr($key1,21,1);
        foreach my $key2 (@protein){
           my $x1 = substr($key2,30,8);
           my $y1 = substr($key2,38,8);
           my $z1 = substr($key2,46,8);
           my $atom2 = substr($key2,13,3);
           my $base = substr($key2,17,3);
           my $numres2 = substr($key2,22,5);
           my $chain2 = substr($key2,21,1);
               $numres2 =~ s/\s+//g;
           my $dis = sqrt(($x - $x1)**2+($y - $y1)**2+($z-$z1)**2);
                        if ($dis <= 4.5 ) {
                            if ($key != $numres2) {
                                 print OUT "$key\t$chain1\t$dis\t$res1\t$atom1\t$numres2\t$chain2\t$base\t$atom2\n";
                            }                     
       }
}
}
close OUT;
}

#########################################################Use hbplus program to calculate hydrogen bond characteristics#################################################################
system "mkdir -p $software_dir/hbplus/job";
system "cp $workdir/output/data_process/$pdb_chain.pdb $software_dir/hbplus/job";
system "cp $workdir/output/data_process/$mut_name.pdb $software_dir/hbplus/job";
chdir "$software_dir/hbplus/job";
system "../hbplus $pdb_chain.pdb";
system "../hbplus $mut_name.pdb";
system "mv ./*.hb2 $workdir/output/bio_feature/HB";
system "rm ./*.pdb";

#####################################################################Calculate NHB features#################################################################################
for my $i (@name){
open(OUT,">","$workdir/output/bio_feature/HB/$i"."_bond.txt") or die "$!";
open(IN,"<","$workdir/output/bio_feature/HB/$i".".hb2") or die "$!";
chomp(my @data=<IN>);
close IN;
foreach my $i(8..$#data) {                
        my $ch1 = substr($data[$i],0,1);      
        my $atom1 = substr($data[$i],1,4);    
        my $res1 = substr($data[$i],6,3);      
        my $ty1 = substr($data[$i],10,3);     
        $atom1 =~ s/^0+//;                    
            
        my $ch2 = substr($data[$i],14,1);     
        my $atom2 = substr($data[$i],15,4);   
        my $res2 = substr($data[$i],20,3);    
        my $ty2 = substr($data[$i],24,3);     
        $atom2 =~ s/^0+//;                    
print OUT "$ch1\t$atom1\t$res1\t$ty1\t$ch2\t$atom2\t$res2\t$ty2\n";          
 }
close OUT;

my @hb;
open(IN,"<", "$workdir/output/bio_feature/HB/$i"."_bond.txt") or die "$!";
while (my $line = <IN>) {
    my @line = split(/\s+/,$line);
    push @hb,"$line[1]_$line[2]";
    push @hb,"$line[5]_$line[6]";   
}
close IN;

my %hash;
foreach my $i (@hb){
    my @num = split(/_/,$i);
if (exists $aa{$num[1]}) {
        $hash{$num[0]}{$num[1]}+=1;
}
}

open(OUT, ">","$workdir/output/bio_feature/HB/$i"."_number.txt") or die "$!";
foreach my $key1 (sort {$a<=>$b} keys %hash){
    my $hash2 = $hash{$key1};
    foreach my $key2 (keys %{$hash2}){
          print OUT "$key1\t$aa{$key2}\t$hash{$key1}{$key2}\n";
    }  
}
close OUT;
}
my @mt_hb;
my @wt_hb;
open(IN,"$workdir/output/bio_feature/HB/$pdb_chain"."_number.txt") or die "$!";
while (my $line = <IN>) {
    my @line = split(/\s+/,$line);
    if ($line[0] == $mutate_pos ) {
        push (@wt_hb,@line)
     }
}
close IN;
if (@wt_hb) {
   
}else{
    push (@wt_hb,0,0,0)  
}
open(INN,"$workdir/output/bio_feature/HB/$mut_name"."_number.txt") or die "$!";
while (my $line = <INN>) {
    my @line = split(/\s+/,$line);
    if ($line[0] == $mutate_pos ) {
        push (@mt_hb,@line)
     }
}
close INN;
if (@mt_hb) {
   
}else{
    push (@mt_hb,0,0,0)  
}
open(OUT,">","$workdir/output/bio_feature/HB/NHB.txt") or die "$!";
foreach my $kk(2..$#wt_hb) {
                my $cha = $mt_hb[$kk] - $wt_hb[$kk];
                print OUT $cha,"\t";
            }
            foreach my $kk(2..$#wt_hb) {
                print OUT $wt_hb[$kk],"\t";}
close OUT;

#########################################################Use naccess program to calculate solvent accessibility characteristics#####################################################
system "mkdir -p $software_dir/naccess/job";
system "cp $workdir/output/data_process/$pdb_chain.pdb $software_dir/naccess/job";
system "cp $workdir/output/data_process/MT_$pdb_chain.pdb $software_dir/naccess/job";
system "cp $workdir/output/data_process/WT_pro.pdb $software_dir/naccess/job";
system "cp $workdir/output/data_process/MT_pro.pdb $software_dir/naccess/job";
chdir "$software_dir/naccess/job";
system "../naccess $pdb_chain.pdb";
system "../naccess MT_$pdb_chain.pdb";
system "../naccess WT_pro.pdb";
system "../naccess MT_pro.pdb";
system "mv $software_dir/naccess/job/MT_$pdb_chain.rsa $workdir/output/bio_feature/ASA";
system "mv $software_dir/naccess/job/MT_pro.rsa $workdir/output/bio_feature/ASA";
system "mv $software_dir/naccess/job/$pdb_chain.rsa $workdir/output/bio_feature/ASA";
system "mv $software_dir/naccess/job/WT_pro.rsa $workdir/output/bio_feature/ASA";
system "rm ./*";

##############################################################################Calculate dASA features################################################################################
my @wt_bind_asa;
my @wt_unbind_asa;
my @mt_bind_asa;
my @mt_unbind_asa;
open(INN,"$workdir/output/bio_feature/ASA/$pdb_chain.rsa") or die "$!";
while (my $line = <INN>) {
    my @line = split(/\s+/,$line);
     if (exists( $aa{$line[1]}) && $line[3]=~ /\d+/ && $line[3] == $mutate_pos) {
        push (@wt_bind_asa,@line[4,6,8,10,12])
     }
}
close INN;
if (@wt_bind_asa) {
   
}else{
    push (@wt_bind_asa,0,0,0,0,0)  
}
open(INN,"$workdir/output/bio_feature/ASA/WT_pro.rsa") or die "$!";
while (my $line = <INN>) {
    my @line = split(/\s+/,$line);
     if (exists( $aa{$line[1]}) && $line[3]=~ /\d+/ && $line[3] == $mutate_pos) {
        push (@wt_unbind_asa,@line[4,6,8,10,12])
     }
}
close INN;
if (@wt_unbind_asa) {
   
}else{
    push (@wt_unbind_asa,0,0,0,0,0)  
}
open(INN,"$workdir/output/bio_feature/ASA/MT_$pdb_chain.rsa") or die "$!";
while (my $line = <INN>) {
    my @line = split(/\s+/,$line);
     if (exists( $aa{$line[1]}) && $line[3]=~ /\d+/ && $line[3] == $mutate_pos) {
        push (@mt_bind_asa,@line[4,6,8,10,12])
     }
}
close INN;
if (@mt_bind_asa) {
   
}else{
    push (@mt_bind_asa,0,0,0,0,0)  
}
open(INN,"$workdir/output/bio_feature/ASA/MT_pro.rsa") or die "$!";
while (my $line = <INN>) {
    my @line = split(/\s+/,$line);
     if (exists( $aa{$line[1]}) && $line[3]=~ /\d+/ && $line[3] == $mutate_pos) {
        push (@mt_unbind_asa,@line[4,6,8,10,12])
     }
}
close INN;
if (@mt_unbind_asa) {
   
}else{
    push (@mt_unbind_asa,0,0,0,0,0)  
}
my @mt_dASA;
my @wt_dASA;
for my $i (0..$#wt_bind_asa){
    my $cha = $wt_unbind_asa[$i] - $wt_bind_asa[$i];
    push (@wt_dASA,$cha)
}
for my $i (0..$#mt_bind_asa){
    my $cha = $mt_unbind_asa[$i] - $mt_bind_asa[$i];
    push (@mt_dASA,$cha)
}

open(OUT,">","$workdir/output/bio_feature/ASA/dASA.txt") or die "$!";
foreach my $kk(0..$#wt_dASA) {
                my $cha = $mt_dASA[$kk] - $wt_dASA[$kk];
                print OUT $cha,"\t";
            }
            foreach my $kk(0..$#wt_dASA) {
                print OUT $wt_dASA[$kk],"\t";}
close OUT;

##############################################################################Calculate IR-dASA features############################################################
open(OUT,">","$workdir/output/bio_feature/ASA/wt_bind_asa.txt") or die "$!";
open(IN,"<","$workdir/output/bio_feature/ASA/$pdb_chain.rsa") or die "$!";
while (my $line = <IN>) {
    my @line = split(/\s+/,$line);
    if ($line[3]=~ /\d+/ && $line[0] =~/RES/ && exists $aa{$line[1]}) {
       print OUT "$line[3]\t$aa{$line[1]}\t$line[4]\t$line[6]\t$line[8]\t$line[10]\t$line[12]\n";
     }
}
close IN;
close OUT;

open(OUT,">","$workdir/output/bio_feature/ASA/wt_unbind_asa.txt") or die "$!";
open(IN,"<","$workdir/output/bio_feature/ASA/WT_pro.rsa") or die "$!";
while (my $line = <IN>) {
    my @line = split(/\s+/,$line);
    if ($line[3]=~ /\d+/ && $line[0] =~/RES/ && exists $aa{$line[1]}) {
       print OUT "$line[3]\t$aa{$line[1]}\t$line[4]\t$line[6]\t$line[8]\t$line[10]\t$line[12]\n";
     }
}
close IN;
close OUT;

open(OUT,">","$workdir/output/bio_feature/ASA/mt_bind_asa.txt") or die "$!";
open(IN,"<","$workdir/output/bio_feature/ASA/MT_$pdb_chain.rsa") or die "$!";
while (my $line = <IN>) {
    my @line = split(/\s+/,$line);
    if ($line[3]=~ /\d+/ && $line[0] =~/RES/ && exists $aa{$line[1]}) {
       print OUT "$line[3]\t$aa{$line[1]}\t$line[4]\t$line[6]\t$line[8]\t$line[10]\t$line[12]\n";
     }
}
close IN;
close OUT;

open(OUT,">","$workdir/output/bio_feature/ASA/mt_unbind_asa.txt") or die "$!";
open(IN,"<","$workdir/output/bio_feature/ASA/MT_pro.rsa") or die "$!";
while (my $line = <IN>) {
    my @line = split(/\s+/,$line);
    if ($line[3]=~ /\d+/ && $line[0] =~/RES/ && exists $aa{$line[1]}) {
       print OUT "$line[3]\t$aa{$line[1]}\t$line[4]\t$line[6]\t$line[8]\t$line[10]\t$line[12]\n";
     }
}
close IN;
close OUT;

open(OUT,">","$workdir/output/bio_feature/ASA/wt_dasa.txt") or die "$!";
open(IN,"<","$workdir/output/bio_feature/ASA/wt_unbind_asa.txt") or die "$!";
chomp(my @wt_pro_asa= <IN>);
close IN;
open(INN,"<","$workdir/output/bio_feature/ASA/wt_bind_asa.txt") or die "$!";
chomp(my @wt_com_asa = <INN>);
close INN;
foreach my $i(0..$#wt_pro_asa) {
            my @aa = split(/\s+/,$wt_pro_asa[$i]);
            my @bb = split(/\s+/,$wt_com_asa[$i]);
            if ($aa[0] eq $bb[0]) {
                print OUT "$aa[0]\t$aa[1]\t";
                foreach my $j(2..$#aa) {
                    my $cha = $aa[$j] - $bb[$j];
                    print OUT "$cha\t";
                }
                print OUT "\n";
            }
}
close OUT;

open(OUT,">","$workdir/output/bio_feature/ASA/mt_dasa.txt") or die "$!";
open(IN,"<","$workdir/output/bio_feature/ASA/mt_unbind_asa.txt") or die "$!";
chomp(my @mt_pro_asa = <IN>);
close IN;
open(INN,"<","$workdir/output/bio_feature/ASA/mt_bind_asa.txt") or die "$!";
chomp(my @mt_com_asa = <INN>);
close INN;
foreach my $i(0..$#mt_pro_asa) {
            my @aa = split(/\s+/,$mt_pro_asa[$i]);
            my @bb = split(/\s+/,$mt_com_asa[$i]);
            if ($aa[0] eq $bb[0]) {
                print OUT "$aa[0]\t$aa[1]\t";
                foreach my $j(2..$#aa) {
                    my $cha = $aa[$j] - $bb[$j];
                    print OUT "$cha\t";
                }
                print OUT "\n";
            }
}
close OUT;

my %inter_pro;
open(INI,"<","$workdir/output/distance/$pdb_chain"."_pro_nul_dis_5.txt") or die "$!";
while (<INI>) {
        chomp;
        my @tem = split(/\s+/);
        if ($tem[2] <= 5) {
            $inter_pro{$tem[0]} = 0;
        }
    }
close INI;

my @mt_inter_dasa = (0,0,0,0,0);
my @wt_inter_dasa = (0,0,0,0,0);

open(INN,"<","$workdir/output/bio_feature/ASA/mt_dasa.txt") or die "$!"; 
while (<INN>) {
        chomp;
        $_ =~ s/^\s+//;
        my @tem = split(/\s+/);
        if (exists $inter_pro{$tem[0]}) {
                foreach my $kk(2..$#tem) {
                    $mt_inter_dasa[$kk-2] += $tem[$kk];
                }
        }
    }
close INN;

open(INN,"<","$workdir/output/bio_feature/ASA/wt_dasa.txt") or die "$!";   
while (<INN>) {
        chomp;
        $_ =~ s/^\s+//;
        my @tem = split(/\s+/);
        if (exists $inter_pro{$tem[0]}) {
                foreach my $kk(2..$#tem) {
                    $wt_inter_dasa[$kk-2] += $tem[$kk];
                }
        }
    }
close INN;

open(OUT,">","$workdir/output/bio_feature/ASA/IR-dASA.txt") or die "$!";
foreach my $kk(0..$#wt_inter_dasa) {
                my $cha = $mt_inter_dasa[$kk] - $wt_inter_dasa[$kk];
                print OUT $cha,"\t";
            }
            foreach my $kk(0..$#wt_inter_dasa) {
                print OUT $wt_inter_dasa[$kk],"\t";}
close OUT;


############################################################################Calculate bRSA features##########################################################
my @wt_com_rsa;
my @mt_com_rsa;
open(INN,"$workdir/output/bio_feature/ASA/$pdb_chain.rsa") or die "$!";
while (my $line = <INN>) {
    my @line = split(/\s+/,$line);
   if (exists( $aa{$line[1]}) && $line[3]=~ /\d+/ && $line[3] == $mutate_pos) {
        push (@wt_com_rsa,@line[5,7,9,11,13])
     }
}
close INN;
if (@wt_com_rsa) {
   
}else{
    push (@wt_com_rsa,0,0,0,0,0)  
}     
open(INN,"$workdir/output/bio_feature/ASA/MT_$pdb_chain.rsa") or die "$!";
while (my $line = <INN>) {
    my @line = split(/\s+/,$line);
    if (exists( $aa{$line[1]}) && $line[3]=~ /\d+/ && $line[3] == $mutate_pos) {
        push (@mt_com_rsa,@line[5,7,9,11,13])
     }
}
close INN;
if (@mt_com_rsa) {
   
}else{
    push (@mt_com_rsa,0,0,0,0,0)  
}   
open(OUT,">","$workdir/output/bio_feature/ASA/bRSA.txt") or die "$!";
foreach my $kk(0..$#wt_com_rsa) {
                my $cha = $mt_com_rsa[$kk] - $wt_com_rsa[$kk];
                print OUT $cha,"\t";
            }
            foreach my $kk(0..$#wt_com_rsa) {
                print OUT $wt_com_rsa[$kk],"\t";}
close OUT;

#############################################################################Calculate uRSA features###########################################################
my @wt_pro_rsa;
my @mt_pro_rsa;
open(INN,"$workdir/output/bio_feature/ASA/WT_pro.rsa") or die "$!";
while (my $line = <INN>) {
    my @line = split(/\s+/,$line);
    if (exists( $aa{$line[1]}) && $line[3]=~ /\d+/ && $line[3] == $mutate_pos) {
        push (@wt_pro_rsa,@line[5,7,9,11,13])
     }
}
close INN;
if (@wt_pro_rsa) {
   
}else{
    push (@wt_pro_rsa,0,0,0,0,0)  
}          
open(INN,"$workdir/output/bio_feature/ASA/MT_pro.rsa") or die "$!";
while (my $line = <INN>) {
    my @line = split(/\s+/,$line);
    if (exists( $aa{$line[1]}) && $line[3]=~ /\d+/ && $line[3] == $mutate_pos) {
        push (@mt_pro_rsa,@line[5,7,9,11,13])
     }
}
close INN;
if (@mt_pro_rsa) {
   
}else{
    push (@mt_pro_rsa,0,0,0,0,0)  
}   
open(OUT,">","$workdir/output/bio_feature/ASA/uRSA.txt") or die "$!";
foreach my $kk(0..$#wt_pro_rsa) {
                my $cha = $mt_pro_rsa[$kk] - $wt_pro_rsa[$kk];
                print OUT $cha,"\t";
            }
            foreach my $kk(0..$#wt_pro_rsa) {
                print OUT $wt_pro_rsa[$kk],"\t";}
close OUT;

################################################################################Calculate JSD features####################################################################################
system"$software_dir/ncbi-blast-2.11.0+/bin/psiblast -query $workdir/output/data_process/WT_pro.fas -db $software_dir/nr/nr  -num_iterations 3  -out 1 -out_ascii_pssm $workdir/output/bio_feature/JSD/pssm/WT.pssm";
system"$software_dir/ncbi-blast-2.11.0+/bin/psiblast -query $workdir//output/data_process/MT_pro.fas -db $software_dir/nr/nr  -num_iterations 3  -out 1 -out_ascii_pssm $workdir/output/bio_feature/JSD/pssm/MT.pssm";
system"/usr/bin/perl $workdir/scripts/4_cscore+_zscore+logical_gzz.pl $workdir";

my %wt_has;
open(IN, "$workdir/output/data_process/WT_pro_seq.txt") or die "$!";
while (my $line = <IN>) {
    chomp ($line);
    my @line = split(/\s+/,$line);
    $wt_has{$line[2]}=$line[0];
    }
close IN;

my $wt_JSD;
open(IN, "$workdir/output/bio_feature/JSD/cscores/WT.txt") or die "$!";
while (my $line = <IN>) {
    chomp ($line);
    if ($line =~/^$wt_has{$mutate_pos}/) {
        my @line = split(/\s+/,$line);
        $wt_JSD=$line[4];
    }
    }
close IN;

my %mt_has;
open(IN, "$workdir/output/data_process/MT_pro_seq.txt") or die "$!";
while (my $line = <IN>) {
    chomp ($line);
    my @line = split(/\s+/,$line);
    $wt_has{$line[2]}=$line[0];
    }
close IN;

my $mt_JSD;
open(IN, "$workdir/output/bio_feature/JSD/cscores/MT.txt") or die "$!";
while (my $line = <IN>) {
    chomp ($line);
    if ($line =~/^$wt_has{$mutate_pos}/) {
        my @line = split(/\s+/,$line);
        $mt_JSD=$line[4];
    }
    }
close IN;

open(OUT, ">$workdir/output/bio_feature/JSD/JSD.txt") or die "$!";
my $cha = $mt_JSD - $wt_JSD;
print OUT "$cha\t$wt_JSD\n";
close OUT;
################################################################################Calculate ENDES features####################################################################################
system "cp $workdir/output/data_process/WT_pro.pdb $software_dir/enrich/unbound";
system "cp $workdir/output/data_process/MT_pro.pdb $software_dir/enrich/unbound";
chdir "$software_dir/enrich";
system "rm energy";
system "./processpdb WT_pro";
system "mv $software_dir/enrich/energy $workdir/output/bio_feature/ENDES/$pdb_chain'.'_enrich.txt";
system "rm energy";
system "./processpdb MT_pro";
system "mv $software_dir/enrich/energy $workdir/output/bio_feature/ENDES/$mut_name'.'_enrich.txt";
system "rm energy";

system "cp $workdir/output/data_process/WT_pro.pdb $software_dir/pinup";
system "cp $workdir/output/data_process/MT_pro.pdb $software_dir/pinup";
chdir "$software_dir/pinup";
system "rm score";
system "rm p001.pdb";
system "mv WT_pro.pdb p001.pdb";
system "./pinup";
system "mv $software_dir/pinup/score $workdir/output/bio_feature/ENDES/$pdb_chain'.'_pinup.txt";
system "rm score";
system "rm p001.pdb";
system "mv MT_pro.pdb p001.pdb";
system "./pinup";
system "mv $software_dir/pinup/score $workdir/output/bio_feature/ENDES/$mut_name'.'_pinup.txt";
system "rm score";
system "rm p001.pdb";

my @wt_endes;
my @mt_endes;
my %wt_endes;
my %wt_endes1;
my %mt_endes;
my %mt_endes1;

open(AA,"<","$workdir/output/bio_feature/ENDES/$pdb_chain"."._pinup.txt") or die "$!";
    while (<AA>) {
        chomp;
        my @temp = split(/\s+/);        
        $wt_endes{$temp[1]} = $temp[3];   
    }
close AA;
if (exists $wt_endes{$mutate_pos}){push @wt_endes,$wt_endes{$mutate_pos}};
open(AA,"<","$workdir/output/bio_feature/ENDES/$pdb_chain"."._enrich.txt") or die "$!"; 
while (<AA>) {
    chomp;
     my @temp = split(/\s+/);          
    $wt_endes1{$temp[1]} = [@temp[3..8]];   
    }
close AA;
if (exists $wt_endes1{$mutate_pos}){push @wt_endes,@{$wt_endes1{$mutate_pos}}};

open(AA,"<","$workdir/output/bio_feature/ENDES/$mut_name"."._pinup.txt") or die "$!";
    while (<AA>) {
        chomp;
        my @temp = split(/\s+/);        
        $mt_endes{$temp[1]} = $temp[3];   
    }
close AA;
if (exists $mt_endes{$mutate_pos}){push @mt_endes,$mt_endes{$mutate_pos}};
open(AA,"<","$workdir/output/bio_feature/ENDES/$mut_name"."._enrich.txt") or die "$!"; 
while (<AA>) {
    chomp;
    my @temp = split(/\s+/);        
    $mt_endes1{$temp[1]} = [@temp[3..8]];   
    }
close AA;
if (exists $mt_endes1{$mutate_pos}){push @mt_endes,@{$mt_endes1{$mutate_pos}}};

open(OUT,">","$workdir/output/bio_feature/ENDES/ENDES.txt") or die "$!";
foreach my $kk(0..$#wt_endes) {
                my $cha = $mt_endes[$kk] - $wt_endes[$kk];
                print OUT $cha,"\t";
            }
            foreach my $kk(0..$#wt_endes) {
                print OUT $wt_endes[$kk],"\t";}
close OUT;
#############################################################################################Connect non-energy features###########################################
my @bio_feature;
open(OUT,">","$workdir/output/bio_feature/bio_feature.txt") or die "$!";
print OUT $pdb_chain."_".$mutate_res;
print OUT "\t";
open(INN,"$workdir/output/bio_feature/ASA/bRSA.txt") or die "$!";
while (my $line = <INN>) {
    chomp $line;
    my @line = split(/\s+/,$line);
    push(@bio_feature,@line)
}
close INN;

open(INN,"$workdir/output/bio_feature/ASA/dASA.txt") or die "$!";
while (my $line = <INN>) {
    chomp $line;
    my @line = split(/\s+/,$line);
    push(@bio_feature,@line)
}
close INN;

open(INN,"$workdir/output/bio_feature/ENDES/ENDES.txt") or die "$!";
while (my $line = <INN>) {
    chomp $line;
    my @line = split(/\s+/,$line);
     push(@bio_feature,@line)
}
close INN;

open(INN,"$workdir/output/output/bio_feature/HB/NHB.txt") or die "$!";
while (my $line = <INN>) {
    chomp $line;
    my @line = split(/\s+/,$line);
     push(@bio_feature,@line)
}
close INN;

open(INN,"$workdir/output/bio_feature/ASA/IR-dASA.txt") or die "$!";
while (my $line = <INN>) {
    chomp $line;
    my @line = split(/\s+/,$line);
     push(@bio_feature,@line)
}
close INN;

open(INN,"$workdir/output/bio_feature/JSD/JSD.txt") or die "$!";
while (my $line = <INN>) {
    chomp $line;
    my @line = split(/\s+/,$line);
     push(@bio_feature,@line)
}
close INN;

open(INN,"$workdir/output/bio_feature/ASA/uRSA.txt") or die "$!";
while (my $line = <INN>) {
    chomp $line;
    my @line = split(/\s+/,$line);
     push(@bio_feature,@line)
}
close INN;

for my $i (@bio_feature){
    print OUT "$i\t";
}

close OUT;

#############################################Use amber software to calculate energy feature files##############################
system "mkdir -p $workdir/output/energy_feature/GB2/wt";
system "cp $workdir/output/data_process/WT_com.pdb $workdir/output/energy_feature/GB2/wt";
system "cp $workdir/output/data_process/WT_dna.pdb $workdir/output/energy_feature/GB2/wt";
system "cp $workdir/output/data_process/WT_pro.pdb $workdir/output/energy_feature/GB2/wt";
system "mv $workdir/output/energy_feature/GB2/wt/WT_com.pdb $workdir/output/energy_feature/GB2/wt/com.pdb";
system "mv $workdir/output/energy_feature/GB2/wt/WT_pro.pdb $workdir/output/energy_feature/GB2/wt/pro.pdb";
system "mv $workdir/output/energy_feature/GB2/wt/WT_dna.pdb $workdir/output/energy_feature/GB2/wt/dna.pdb";
system "cp $workdir/scripts/GB2/* $workdir/output/energy_feature/GB2/wt";
chdir "$workdir/output/energy_feature/GB2/wt";
system "sh $workdir/output/energy_feature/GB2/wt/amber.sh";
system "mv $workdir/output/energy_feature/GB2/wt/Final_result_mmpbsa_GB2_pair.dat $workdir/output/energy_feature/GB2/$pdb_chain'.'_pair.dat";

system "mkdir -p $workdir/output/energy_feature/GB2/mt";
system "cp $workdir/output/data_process/COM.pdb $workdir/output/energy_feature/GB2/mt";
system "cp $workdir/output/data_process/DNA.pdb $workdir/output/energy_feature/GB2/mt";
system "cp $workdir/output/data_process/MT_pro.pdb /$workdir/output/energy_feature/GB2/mt";
system "mv $workdir/output/energy_feature/GB2/mt/COM.pdb $workdir/output/energy_feature/GB2/mt/com.pdb";
system "mv $workdir/output/energy_feature/GB2/mt/MT_pro.pdb $workdir/output/energy_feature/GB2/mt/pro.pdb";
system "mv $workdir/output/energy_feature/GB2/mt/DNA.pdb $workdir/output/energy_feature/GB2/mt/dna.pdb";
system "cp $workdir/scripts/GB2/* $workdir/output/energy_feature/GB2/mt";
chdir "$workdir/output/energy_feature/GB2/mt";
system "sh $workdir/output/energy_feature/GB2/mt/amber.sh";
system "mv $workdir/output/energy_feature/GB2/mt/Final_result_mmpbsa_GB2_pair.dat $workdir/output/energy_feature/GB2/$mut_name'.'_pair.dat";

open(IN,"<","$workdir/output/energy_feature/GB2/$pdb_chain"."._pair.dat") or die "$!";
my @wt_data = <IN>;
close IN;

my @wt;
foreach my $j(81..84) {       
    my @tem = split(/\s+/,$wt_data[$j]);
    push @wt,$tem[1];
        }
my @wt_aaa = split(/\s+/,$wt_data[89]);  
push @wt,$wt_aaa[2];

open(IN,"<","$workdir/output/energy_feature/GB2/$mut_name"."._pair.dat") or die "$!";
my @mt_data = <IN>;
close IN;

my @mt;
foreach my $j(81..84) {       
    my @tem = split(/\s+/,$mt_data[$j]);
    push @mt,$tem[1];
        }
my @mt_aaa = split(/\s+/,$mt_data[89]);  
push @mt,$mt_aaa[2];

open(OUT,">","$workdir/output/energy_feature/GB2/EWC.txt") or die "$!";
foreach my $kk(0..$#wt) {
                my $cha = $mt[$kk] - $wt[$kk];
                print OUT $cha,"\t";
            }
            foreach my $kk(0..$#wt) {
                print OUT $wt[$kk],"\t";}
close OUT;

####result.txt####
open(OUT,">","$workdir/output/result.txt") or die "$!";
print OUT "PDB ID:$pdb_name\nChian ID:$chain_name\nPosition:$mutate_res";
close OUT;
########################Random forest prediction#############################
system "/home/yaojiang/anaconda3/bin/python $workdir/scripts/MPR-regression.py $workdir";

########################Random forest classification#############################
system "/home/yaojiang/anaconda3/bin/python $workdir/scripts/MPR-classification.py $workdir";

