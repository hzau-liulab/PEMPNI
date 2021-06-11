use warnings;
use strict;

#Examples:perl protein_DNA.pl 1AAY A D120A
#Notice:You need to set the correct working path and software installation path.
#And you need to place the wild structure file in the ./WT_PDB directory.

#############################################################Set working directory#################################################
my $workdir="/home/yaojiang/protein-DNA/";                #Set to your own working directory 
my $software_dir="/var/www/html/PEMPNI/software/";        #Set the software installation directory

my %aa=('ALA'=>'A','ARG'=>'R','ASN'=>'N','ASP'=>'D','CYS'=>'C','GLN'=>'Q','GLU'=>'E','GLY'=>'G','HIS'=>'H','ILE'=>'I',
        'LEU'=>'L','LYS'=>'K','MET'=>'M','PHE'=>'F','PRO'=>'P','SER'=>'S','THR'=>'T','TRP'=>'W','TYR'=>'Y','VAL'=>'V');
my %reverse_aa=reverse %aa;

my $pdb_name = $ARGV[0];                       #1AAY
my $chain_name = $ARGV[1];                     #A
my $mutate_res = $ARGV[2];                     #D120A
my $mutated_res=substr($mutate_res,-1,1);
my $mutate_restyp = $reverse_aa{$mutated_res}; #ALA
my $pdb_chain = $pdb_name.'_'.$chain_name;     #1AAY_A
my $mut_name="$pdb_chain.mut.$mutate_res";     #1AAY_A.mut.D120A
my $mutate_pos =$mutate_res;
   $mutate_pos =~s/[A-Z]//g;                   #120
my $wild_res=substr($mutate_res,0,1);          #D
my $wlid_restyp=$reverse_aa{$wild_res};        #ASP

my @name;
push(@name,$pdb_chain);
push(@name,$mut_name);

system "mkdir -p $workdir/$mut_name/data_process";
system "mkdir -p $workdir/$mut_name/energy_feature";
system "mkdir -p $workdir/$mut_name/energy_feature/GB1";
system "mkdir -p $workdir/$mut_name/distance";
system "mkdir -p $workdir/$mut_name/bio_feature";
system "mkdir -p $workdir/$mut_name/bio_feature/ASA";
system "mkdir -p $workdir/$mut_name/bio_feature/contact";
system "mkdir -p $workdir/$mut_name/bio_feature/HB";
system "cp $workdir/WT_PDB/$pdb_name.pdb $workdir/$mut_name/data_process";

############################################################Original structure processing#############################################
open(IN,"$workdir/$mut_name/data_process/$pdb_name.pdb") or die "can not open this files $!\n"; 
open(OUT,">$workdir/$mut_name/data_process/WT_pro.pdb") or die "can not open $!\n";
open(OUTTTTT,">$workdir/$mut_name/data_process/WT_dna.pdb") or die "can not open $!\n";
open(OUTT,">$workdir/$mut_name/data_process/pre_dna.pdb") or die "can not open $!\n";
open(OUTTT,">$workdir/$mut_name/data_process/WT_com.pdb") or die "can not open $!\n";
open(OUTTTT,">$workdir/$mut_name/data_process/$pdb_chain.pdb") or die "can not open $!\n";
while (my $line=<IN>) {
    my $atom = substr($line,77,1);
    $atom =~s/\s+//g;
    if($line=~/^ATOM/){  
    if ($line=~/^(ATOM).*D[AGCT]/) {print OUTT "$line"; }
    elsif( $ARGV[1] eq substr($line,21,1) and $atom ne 'H') {print OUT "$line";}
     }
    }
close IN;
close OUT;
close OUTT;
my %hhh;
open(IN, "$workdir/$mut_name/data_process/WT_pro.pdb") or die "can not open $!\n";
my @propro=<IN>;
close IN;

open(IN, "$workdir/$mut_name/data_process/pre_dna.pdb") or die "can not open $!\n";
my @nulnul=<IN>;
close IN;

open(OOUT,">$workdir/$mut_name/data_process/pree_dna.pdb") or die "can not open $!\n";

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
                        if ($dis <= 5) {$hhh{$chain1}+=1;}
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

open(INN, "$workdir/$mut_name/data_process/pree_dna.pdb") or die "can not open $!\n";
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

open(IN, "$workdir/$mut_name/data_process/WT_dna.pdb") or die "can not open $!\n";
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

############################################Check whether the file meets the amber calculation requirements#####################
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
my $out = "$workdir/$mut_name/data_process/";
system "$software_dir/modeller9.25/bin/mod9.25 $workdir/scripts/mutate_model.py $out $mutate_res $pdb_chain $mutate_pos $mutate_restyp $chain_name >>./log.txt";

open(IN,"$workdir/$mut_name/data_process/$mut_name.pdb") or die "can not open this files\n"; 
open(OUT,">$workdir/$mut_name/data_process/MT_pro.pdb") or die "can not open $!\n";
open(OUTT,">$workdir/$mut_name/data_process/MT_$pdb_chain.pdb") or die "can not open $!\n";
while (my $line=<IN>) {
    if($line=~/^ATOM/){  
    if ($line=~/^(ATOM).*D[AGCT]/) {print OUTT "$line";}
    elsif( $ARGV[1] eq substr($line,21,1)) {print OUT "$line";print OUTT "$line";}
     }
    }
close IN;
close OUT;
close OUTT;

open(OUTTTTT,">$workdir/$mut_name/data_process/MT_dna.pdb") or die "can not open $!\n";
my @nul;

open(INN,"$workdir/$mut_name/data_process/$mut_name.pdb") or die "can not open this files\n"; 
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
    push(@nul,$pdb[$i]);
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

open(OU,">$workdir/$mut_name/data_process/MT_com.pdb") or die "can not open $!\n";
open(IN,"$workdir/$mut_name/data_process/MT_dna.pdb") or die "can not open $!\n";
while (my $line = <IN>) {
    chomp $line;
    print OU "$line\n";
}
close IN;

open(IN,"$workdir/$mut_name/data_process/MT_pro.pdb") or die "can not open $!\n";
while (my $line = <IN>) {
    chomp $line;
    print OU "$line\n";
}
close IN;
close OU;

##############################################################calculate diff—distance files#########################################################
my $tem = $pdb_name."_".$chain_name."_diff-distance";  
my @te = split(/_/,$tem);                 
my $site = $mutate_pos;                   #120              
my $wt = substr($mutate_res,0,1);         #D
my @pro_coor;                             
my %pdb_coor;                           
open(INN,"<","$workdir/$mut_name/data_process/WT_com.pdb") or die "$!";    #1AAY_A.pdb 
while (<INN>) {
    chomp;
    if ($_ =~ /^ATOM/) {
        my $num = substr($_,22,5);
           $num =~ s/\s+//g;
        my $res = substr($_,17,3);
        my $cha = substr($_,21,1);
        if ($num eq $site and $cha eq $te[1] and $aa{$res} eq $wt) {      
                push @pro_coor,$_;                                       
        }else{
            if (exists $pdb_coor{"$num\t$res\t$cha"}) {              
                push @{$pdb_coor{"$num\t$res\t$cha"}},$_;
                }else{
                    my @bb;
                    push @bb,$_;
                    $pdb_coor{"$num\t$res\t$cha"} = \@bb;
                    }
                }
            }
        }
        close INN;
        
        my %inter_type;   
        foreach my $kk(sort keys%pdb_coor) {
            my @distance;
            foreach my $k1(@{$pdb_coor{$kk}}) {          
                my $x = substr($k1,30,8);
                my $y = substr($k1,38,8);
                my $z = substr($k1,46,8);     
                foreach my $k2(@pro_coor) {
                    my $x1 = substr($k2,30,8);
                    my $y1 = substr($k2,38,8);
                    my $z1 = substr($k2,46,8);
                    my $disan = sqrt(($x - $x1)**2+($y - $y1)**2+($z-$z1)**2);
                    push @distance,$disan;
                }
            }  
            @distance = sort{$a <=> $b}@distance;
            my $dis = $distance[0];              
            if ($dis <= 3) {
                $inter_type{$kk} = "$dis\tH-bond";      
            }elsif($dis > 3 and $dis <= 3.5) {
                $inter_type{$kk} = "$dis\tstacking";
            }elsif($dis > 3.5 and $dis <= 4) {
                $inter_type{$kk} = "$dis\twater-mediated";
            }elsif($dis > 4 and $dis <= 5) {
                $inter_type{$kk} = "$dis\thydrophobic";
            }else{
                $inter_type{$kk} = "$dis\tother";
            }        
        }      
open(OUT,">","$workdir/$mut_name/distance/$tem.txt") or die "$!";   
foreach my $key(sort keys%inter_type) {
        my @key=split(/\s+/,$key);
        if (exists $aa{$key[1]}) {
        print OUT "$key\t$inter_type{$key}\n";
    }
}
close OUT;

#################################Calculate the file whose protein-nucleic acid atom distance is less than 5#################################################
for my $i (@name){
my @nucleic_acid;
my @protein;

my @amino_acid = keys %aa;
open(IN,"<","$workdir/$mut_name/data_process/$i.pdb") or die "$!";
while (my $line = <IN>) {
    if ($line =~/^ATOM/) {
    my $aa=substr($line,17,3);
    $aa =~s/\s+//g; 
        if (grep {$_ eq $aa}@amino_acid) {
            push (@protein,$line);
        }else{
            push (@nucleic_acid,$line);
            }
    }else{
        next;
    }     
}
close IN;

open(OUT,">","$workdir/$mut_name/distance/$i"."_pro_nul_dis_5.txt") or die "$!";
foreach my $key1 (@protein){
    my $numres1 = substr($key1,22,5);
    $numres1=~s/\s+//g;
    my $x = substr($key1,30,8);
    my $y = substr($key1,38,8);
    my $z = substr($key1,46,8);
    my $atom1 = substr($key1,13,3);
    my $res1 = $aa{substr($key1,17,3)};
    my $chain1 = substr($key1,21,1);
        foreach my $key2 (@nucleic_acid){
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
                            print OUT "$numres1\t$chain1\t$dis\t$res1\t$atom1\t$numres2\t$chain2\t$base\t$atom2\n";
                        }else{next;}
        }
    }
close OUT;
#################################Calculate files with protein-protein atom distance less than 4.5#################################################
open(OUT,">","$workdir/$mut_name/distance/$i"."_pro_pro_dis_4.5.txt") or die "$!";
foreach my $key1 (@protein){
    my $numres1 = substr($key1,22,5);
    $numres1=~s/\s+//g;
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
                            if ($numres1 ne $numres2) {
                                 print OUT "$numres1\t$chain1\t$dis\t$res1\t$atom1\t$numres2\t$chain2\t$base\t$atom2\n";
                            }                     
       }
}
}
close OUT;
}

####################################################################Calculate CFNA features####################################################
for my $i (@name){
my %res_atom;
my %res_contact;
open(IN,"<","$workdir/$mut_name/distance/$i"."_pro_nul_dis_5.txt") or die "$!";
while (<IN>) {
    chomp;
    my @tem = split(/\s+/);
    if ($tem[2] < 4.5) {                              
        $res_atom{$tem[0]} ++;
        $res_contact{$tem[0]}{$tem[5]} ++;
        }
    }
close IN;
open(OUT,">","$workdir/$mut_name/bio_feature/contact/$i"."_contact_nul.txt") or die "$!";
foreach my $key(sort{$a <=> $b}keys%res_atom) {
    print OUT "$key\t";
    my $contact = scalar keys%{$res_contact{$key}};  
    my $density = $res_atom{$key}/$contact;          
    print OUT "$contact\t$density\n";
}
close OUT;
}

my @mt_contact_nul;
my @wt_contact_nul;
open(IN,"$workdir/$mut_name/bio_feature/contact/$pdb_chain"."_contact_nul.txt") or die "$!";
while (my $line = <IN>) {
    my @line = split(/\s+/,$line);
    if ($line[0] eq $mutate_pos ) {
        push (@wt_contact_nul,@line);
     }
}
close IN;

if (@wt_contact_nul) {
   
}else{
    push (@wt_contact_nul,0,0,0);  
}

open(INN,"$workdir/$mut_name/bio_feature/contact/$mut_name"."_contact_nul.txt") or die "$!";
while (my $line = <INN>) {
    my @line = split(/\s+/,$line);
    if ($line[0] eq $mutate_pos ) {
        push (@mt_contact_nul,@line);
     }
}
close INN;
if (@mt_contact_nul) {
   
}else{
    push (@mt_contact_nul,0,0,0);  
}
open(OUT,">","$workdir/$mut_name/bio_feature/contact/CFNA.txt") or die "$!";
    
foreach my $kk(1..$#wt_contact_nul) {
                my $cha = $mt_contact_nul[$kk] - $wt_contact_nul[$kk];
                print OUT $cha,"\t";
            }
            foreach my $kk(1..$#wt_contact_nul) {
                print OUT $wt_contact_nul[$kk],"\t";}
close OUT;

###################################################################Calculate CFAA features####################################################
for my $i (@name){
my %res_atom;
my %res_contact;
open(IN,"<","$workdir/$mut_name/distance/$i"."_pro_pro_dis_4.5.txt") or die "$!";
while (<IN>) {
    chomp;
    my @tem = split(/\s+/);
    if ($tem[2] < 4.5) {                              
    $res_atom{$tem[0]} ++;
    $res_contact{$tem[0]}{$tem[5]} ++;
    }
}
close IN;
open(OUT,">","$workdir/$mut_name/bio_feature/contact/$i"."_contact_pro.txt") or die "$!";
foreach my $key(sort{$a <=> $b}keys%res_atom) {
    print OUT "$key\t";
    my $contact = scalar keys%{$res_contact{$key}};  
    my $density = $res_atom{$key}/$contact;          
    print OUT "$contact\t$density\n";
    }
close OUT;
}

my @mt_contact_pro;
my @wt_contact_pro;
open(IN,"$workdir/$mut_name/bio_feature/contact/$pdb_chain"."_contact_pro.txt") or die "$!";
while (my $line = <IN>) {
    my @line = split(/\s+/,$line);
    if ($line[0] eq $mutate_pos ) {
        push (@wt_contact_pro,@line);
     }
}
close IN;
if (@wt_contact_pro) {
   
}else{
    push (@wt_contact_pro,0,0,0);  
}
open(INN,"$workdir/$mut_name/bio_feature/contact/$mut_name"."_contact_pro.txt") or die "$!";
while (my $line = <INN>) {
    my @line = split(/\s+/,$line);
    if ($line[0] eq $mutate_pos ) {
        push (@mt_contact_pro,@line);
     }
}
close INN;
if (@mt_contact_pro) {
   
}else{
    push (@mt_contact_pro,0,0,0);  
}
open(OUT,">","$workdir/$mut_name/bio_feature/contact/CFAA.txt") or die "$!";
foreach my $kk(1..$#wt_contact_pro) {
                my $cha = $mt_contact_pro[$kk] - $wt_contact_pro[$kk];
                print OUT $cha,"\t";
            }
            foreach my $kk(1..$#wt_contact_pro) {
                print OUT $wt_contact_pro[$kk],"\t";}
close OUT;

#################################################################Calculate IR-CFAA features#########################################
my %inter_pro;
open(INI,"<","$workdir/$mut_name/distance/$pdb_chain"."_pro_nul_dis_5.txt") or die "$!";   
while (<INI>) {
    chomp;
    my @tem = split(/\s+/);
    if ($tem[2] <= 5) {
        $inter_pro{$tem[0]} = 0;
        }
    }
    close INI;
    
my @inter_contact_pro_wt = (0,0);
my @inter_contact_pro_mt = (0,0);
open(INN,"<","$workdir/$mut_name/bio_feature/contact/$pdb_chain"."_contact_pro.txt") or die "$!";
 while (<INN>) {
        chomp;
        $_ =~ s/^\s+//;
        my @tem = split(/\s+/);
        if (exists $inter_pro{$tem[0]}) {
                foreach my $kk(1..$#tem) {
                    $inter_contact_pro_wt[$kk-1] += $tem[$kk];
                }     
        }
    }
close INN;
open(INN,"<","$workdir/$mut_name/bio_feature/contact/$mut_name"."_contact_pro.txt") or die "$!";
while (<INN>) {
        chomp;
        $_ =~ s/^\s+//;
        my @tem = split(/\s+/);
        if (exists $inter_pro{$tem[0]}) {
                foreach my $kk(1..$#tem) {
                    $inter_contact_pro_mt[$kk-1] += $tem[$kk];
                }     
        }
    }
close INN;
open(OUT,">","$workdir/$mut_name/bio_feature/contact/IR-CFAA.txt") or die "$!";
foreach my $kk(0..$#inter_contact_pro_wt) {
                my $cha = $inter_contact_pro_mt[$kk] - $inter_contact_pro_wt[$kk];
                print OUT $cha,"\t";
            }
            foreach my $kk(0..$#inter_contact_pro_wt) {
                print OUT $inter_contact_pro_wt[$kk],"\t";}
close OUT;

#########################################################Use hbplus program to calculate hydrogen bond characteristics####################################################
system "mkdir -p $software_dir/hbplus/job";
system "cp $workdir/$mut_name/data_process/$pdb_chain.pdb $software_dir/hbplus/job";
system "cp $workdir/$mut_name/data_process/$mut_name.pdb $software_dir/hbplus/job";
chdir "$software_dir/hbplus/job";
system "../hbplus $pdb_chain.pdb";
system "../hbplus $mut_name.pdb";
system "mv ./*.hb2 $workdir/$mut_name/bio_feature/HB";
system "rm ./*.pdb";

###############################################################Calculate NHB features#####################################################################################
for my $i (@name){
open(OUT,">","$workdir/$mut_name/bio_feature/HB/$i"."_bond.txt") or die "$!";
open(IN,"<","$workdir/$mut_name/bio_feature/HB/$i".".hb2") or die "$!";
chomp(my @data=<IN>);
close IN;
foreach my $i(8..$#data) {                #A0170-ARG NH2 B0002- DG OP2 2.55 SH  -2 -1.00 112.7  1.99 106.6 117.4     1
        my $ch1 = substr($data[$i],0,1);      #A
        my $atom1 = substr($data[$i],1,4);    #0170
        my $res1 = substr($data[$i],6,3);     #ARG  
        my $ty1 = substr($data[$i],10,3);     #NH2
        $atom1 =~ s/^0+//;                    #170
            
        my $ch2 = substr($data[$i],14,1);     #B
        my $atom2 = substr($data[$i],15,4);   #00002
        my $res2 = substr($data[$i],20,3);    #DG
        my $ty2 = substr($data[$i],24,3);     #OP2
        $atom2 =~ s/^0+//;                    #2
print OUT "$ch1\t$atom1\t$res1\t$ty1\t$ch2\t$atom2\t$res2\t$ty2\n";  #A	170	ARG	NH2	B	2	 DG	OP2            
 }
close OUT;

my @hb;
open(IN,"<", "$workdir/$mut_name/bio_feature/HB/$i"."_bond.txt") or die "$!";
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

open(OUT, ">","$workdir/$mut_name/bio_feature/HB/$i"."_number.txt") or die "$!";
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
open(IN,"$workdir/$mut_name/bio_feature/HB/$pdb_chain"."_number.txt") or die "$!";
while (my $line = <IN>) {
    my @line = split(/\s+/,$line);
    if ($line[0] eq $mutate_pos ) {
        push (@wt_hb,@line);
     }
}
close IN;
if (@wt_hb) {
   
}else{
    push (@wt_hb,0,0,0);  
}
open(INN,"$workdir/$mut_name/bio_feature/HB/$mut_name"."_number.txt") or die "$!";
while (my $line = <INN>) {
    my @line = split(/\s+/,$line);
    if ($line[0] eq $mutate_pos ) {
        push (@mt_hb,@line);
     }
}
close INN;
if (@mt_hb) {
   
}else{
    push (@mt_hb,0,0,0);  
}
open(OUT,">","$workdir/$mut_name/bio_feature/HB/NHB.txt") or die "$!";
foreach my $kk(2..$#wt_hb) {
                my $cha = $mt_hb[$kk] - $wt_hb[$kk];
                print OUT $cha,"\t";
            }
            foreach my $kk(2..$#wt_hb) {
                print OUT $wt_hb[$kk],"\t";}
close OUT;

########################################Use naccess program to calculate solvent accessibility characteristics#################################################
system "mkdir -p $software_dir/naccess/job";
system "cp $workdir/$mut_name/data_process/$pdb_chain.pdb $software_dir/naccess/job";
system "cp $workdir/$mut_name/data_process/MT_$pdb_chain.pdb $software_dir/naccess/job";
system "cp $workdir/$mut_name/data_process/WT_pro.pdb $software_dir/naccess/job";
system "cp $workdir/$mut_name/data_process/MT_pro.pdb $software_dir/naccess/job";
chdir "$software_dir/naccess/job";
system "../naccess $pdb_chain.pdb";
system "../naccess MT_$pdb_chain.pdb";
system "../naccess WT_pro.pdb";
system "../naccess MT_pro.pdb";
system "mv $software_dir/naccess/job/MT_$pdb_chain.rsa $workdir/$mut_name/bio_feature/ASA";
system "mv $software_dir/naccess/job/MT_pro.rsa $workdir/$mut_name/bio_feature/ASA";
system "mv $software_dir/naccess/job/$pdb_chain.rsa $workdir/$mut_name/bio_feature/ASA";
system "mv $software_dir/naccess/job/WT_pro.rsa $workdir/$mut_name/bio_feature/ASA";
system "rm ./*";

#######################################################Calculate dASA features#######################################################
my @wt_bind_asa;
my @wt_unbind_asa;
my @mt_bind_asa;
my @mt_unbind_asa;
open(INN,"$workdir/$mut_name/bio_feature/ASA/$pdb_chain.rsa") or die "$!";
while (my $line = <INN>) {
    my @line = split(/\s+/,$line);
     if (exists( $aa{$line[1]}) && $line[3]=~ /\d+/ && $line[3] eq $mutate_pos) {
        push (@wt_bind_asa,@line[4,6,8,10,12]);
     }
}
close INN;
if (@wt_bind_asa) {
   
}else{
    push (@wt_bind_asa,0,0,0,0,0);  
}
open(INN,"$workdir/$mut_name/bio_feature/ASA/WT_pro.rsa") or die "$!";
while (my $line = <INN>) {
    my @line = split(/\s+/,$line);
     if (exists( $aa{$line[1]}) && $line[3]=~ /\d+/ && $line[3] eq $mutate_pos) {
        push (@wt_unbind_asa,@line[4,6,8,10,12]);
     }
}
close INN;
if (@wt_unbind_asa) {
   
}else{
    push (@wt_unbind_asa,0,0,0,0,0);  
}

open(INN,"$workdir/$mut_name/bio_feature/ASA/MT_$pdb_chain.rsa") or die "$!";
while (my $line = <INN>) {
    my @line = split(/\s+/,$line);
     if (exists( $aa{$line[1]}) && $line[3]=~ /\d+/ && $line[3] eq $mutate_pos) {
        push (@mt_bind_asa,@line[4,6,8,10,12]);
     }
}
close INN;
if (@mt_bind_asa) {
   
}else{
    push (@mt_bind_asa,0,0,0,0,0);  
}
open(INN,"$workdir/$mut_name/bio_feature/ASA/MT_pro.rsa") or die "$!";
while (my $line = <INN>) {
    my @line = split(/\s+/,$line);
     if (exists( $aa{$line[1]}) && $line[3]=~ /\d+/ && $line[3] eq $mutate_pos) {
        push (@mt_unbind_asa,@line[4,6,8,10,12]);
     }
}
close INN;
if (@mt_unbind_asa) {
   
}else{
    push (@mt_unbind_asa,0,0,0,0,0);  
}

my @mt_dASA;
my @wt_dASA;
for my $i (0..$#wt_bind_asa){
    my $cha = $wt_unbind_asa[$i] - $wt_bind_asa[$i];
    push (@wt_dASA,$cha);
}
for my $i (0..$#mt_bind_asa){
    my $cha = $mt_unbind_asa[$i] - $mt_bind_asa[$i];
    push (@mt_dASA,$cha);
}

open(OUT,">","$workdir/$mut_name/bio_feature/ASA/dASA.txt") or die "$!";
foreach my $kk(0..$#wt_dASA) {
                my $cha = $mt_dASA[$kk] - $wt_dASA[$kk];
                print OUT $cha,"\t";
            }
            foreach my $kk(0..$#wt_dASA) {
                print OUT $wt_dASA[$kk],"\t";}
close OUT;

###########################################################Connect non-energy features############################################################
my @bio_feature;
open(OUT,">","$workdir/$mut_name/bio_feature/bio_feature.txt") or die "$!";
print OUT $pdb_chain."_".$mutate_res;
print OUT "\t";
open(INN,"$workdir/$mut_name/bio_feature/contact/CFNA.txt") or die "$!";
while (my $line = <INN>) {
    chomp $line;
    my @line = split(/\s+/,$line);
    push(@bio_feature,@line);
}
close INN;

open(INN,"$workdir/$mut_name/bio_feature/contact/CFAA.txt") or die "$!";
while (my $line = <INN>) {
    chomp $line;
    my @line = split(/\s+/,$line);
    push(@bio_feature,@line);
}
close INN;

open(INN,"$workdir/$mut_name/bio_feature/ASA/dASA.txt") or die "$!";
while (my $line = <INN>) {
    chomp $line;
    my @line = split(/\s+/,$line);
     push(@bio_feature,@line);
}
close INN;

open(INN,"$workdir/$mut_name/bio_feature/HB/NHB.txt") or die "$!";
while (my $line = <INN>) {
    chomp $line;
    my @line = split(/\s+/,$line);
     push(@bio_feature,@line);
}
close INN;

open(INN,"$workdir/$mut_name/bio_feature/contact/IR-CFAA.txt") or die "$!";
while (my $line = <INN>) {
    chomp $line;
    my @line = split(/\s+/,$line);
     push(@bio_feature,@line);
}
close INN;

for my $i (@bio_feature){
    print OUT "$i\t";
}

close OUT;

##############################################Use amber software to calculate energy feature files##############################
system "mkdir -p $workdir/$mut_name/energy_feature/GB1/wt";
system "cp $workdir/$mut_name/data_process/WT_com.pdb $workdir/$mut_name/energy_feature/GB1/wt";
system "cp $workdir/$mut_name/data_process/WT_dna.pdb $workdir/$mut_name/energy_feature/GB1/wt";
system "cp $workdir/$mut_name/data_process/WT_pro.pdb $workdir/$mut_name/energy_feature/GB1/wt";
system "mv $workdir/$mut_name/energy_feature/GB1/wt/WT_com.pdb $workdir/$mut_name/energy_feature/GB1/wt/com.pdb";
system "mv $workdir/$mut_name/energy_feature/GB1/wt/WT_pro.pdb $workdir/$mut_name/energy_feature/GB1/wt/pro.pdb";
system "mv $workdir/$mut_name/energy_feature/GB1/wt/WT_dna.pdb $workdir/$mut_name/energy_feature/GB1/wt/dna.pdb";
system "cp $workdir/scripts/GB1/* $workdir/$mut_name/energy_feature/GB1/wt";
chdir "$workdir/$mut_name/energy_feature/GB1/wt";
system "sh $workdir/$mut_name/energy_feature/GB1/wt/amber.sh";
system "mv $workdir/$mut_name/energy_feature/GB1/wt/FINAL_DECOMP_MMPBSA_GB1_pair.dat $workdir/$mut_name/energy_feature/GB1/$pdb_chain'.'_pair.dat";

system "mkdir -p $workdir/$mut_name/energy_feature/GB1/mt";
system "cp $workdir/$mut_name/data_process/MT_com.pdb $workdir/$mut_name/energy_feature/GB1/mt";
system "cp $workdir/$mut_name/data_process/MT_dna.pdb $workdir/$mut_name/energy_feature/GB1/mt";
system "cp $workdir/$mut_name/data_process/MT_pro.pdb $workdir/$mut_name/energy_feature/GB1/mt";
system "mv $workdir/$mut_name/energy_feature/GB1/mt/MT_com.pdb $workdir/$mut_name/energy_feature/GB1/mt/com.pdb";
system "mv $workdir/$mut_name/energy_feature/GB1/mt/MT_pro.pdb $workdir/$mut_name/energy_feature/GB1/mt/pro.pdb";
system "mv $workdir/$mut_name/energy_feature/GB1/mt/MT_dna.pdb $workdir/$mut_name/energy_feature/GB1/mt/dna.pdb";
system "cp $workdir/scripts/GB1/* $workdir/$mut_name/energy_feature/GB1/mt";
chdir "$workdir/$mut_name/energy_feature/GB1/mt";
system "sh $workdir/$mut_name/energy_feature/GB1/mt/amber.sh";
system "mv $workdir/$mut_name/energy_feature/GB1/mt/FINAL_DECOMP_MMPBSA_GB1_pair.dat $workdir/$mut_name/energy_feature/GB1/$mut_name'.'_pair.dat";

system "cp $workdir/$mut_name/data_process/MT_com.pdb $workdir/$mut_name/energy_feature/GB1/";
system "mv $workdir/$mut_name/energy_feature/GB1/MT_com.pdb $workdir/$mut_name/energy_feature/GB1/$mut_name'.'_com.pdb";
system "cp $workdir/$mut_name/data_process/WT_com.pdb $workdir/$mut_name/energy_feature/GB1/";
system "mv $workdir/$mut_name/energy_feature/GB1/WT_com.pdb $workdir/$mut_name/energy_feature/GB1/$pdb_chain'.'_com.pdb";

###################################Process the residue pair energy decomposition file output by amber#########################
for my $k (@name){
my @data;
open(IN,"<","$workdir/$mut_name/energy_feature/GB1/$k"."._pair.dat") or die "$!";
while (<IN>) {
    chomp;
    if ($_ =~ /S,i,d,e,c,h,a,i,n/) {last;
    }else{
       push @data,$_; 
      }
    }
close IN;
            
my %set;
foreach my $i(9..($#data-1)) {
    my @tem = split(/,/,$data[$i]);
    $set{$tem[0]} = $i;
    }
my @set_index = sort{$set{$a} <=> $set{$b}} keys %set;
            
my %res;
my $m = 0;
open(INN,"<","$workdir/$mut_name/energy_feature/GB1/$k"."._com.pdb") or die "$!";
while (<INN>) {
    chomp;
if ($_ =~ /^ATOM/) {
    my $residue = substr($_,17,3);
    my $chain = substr($_,21,1);
    my $number = substr($_,22,5);
    $m ++;
    $res{"$residue\t$chain\t$number"} = $m;
    }
}
close INN;
my @res_index = sort{$res{$a} <=> $res{$b}} keys %res;
            
my %map;
foreach my $k(0..$#set_index) {
     $map{$set_index[$k]} = $res_index[$k];
}
open(OUT,">","$workdir/$mut_name/energy_feature/GB1/$k"."_pair_out.txt") or die "$!";
foreach my $k(9..($#data-1)) {
my @temp = split(/,/,$data[$k]);
print OUT "$map{$temp[0]}\t$map{$temp[1]}\t";
print OUT "$temp[2]\t$temp[5]\t$temp[8]\t$temp[11]\t$temp[14]\t$temp[17]\n";
}
close OUT;
}

############################################Calculate EPI features#######################################
open(OUT,">","$workdir/$mut_name/energy_feature/GB1/GB1_WT_EPI.txt") or die "$!"; 
my %wt_inter_pro;            
my %wt_inter_nul;            
open(INI,"<","$workdir/$mut_name/distance/$pdb_chain"."_pro_nul_dis_5.txt") or die "$!";
while (<INI>) {
        chomp;
        my @tem = split(/\s+/);
        if ($tem[2] <= 5) {
            $wt_inter_pro{$tem[0]} = 0;
            $wt_inter_nul{"$tem[6]\t$tem[5]"} = 0;
        }
    }
    close INI;
    
    my @wt_inter_energy = (0,0,0,0,0);    
    my @wt_pro_energy = (0,0,0,0,0);      
    my @wt_nul_energy = (0,0,0,0,0);      
    
open(INN,"<","$workdir/$mut_name/energy_feature/GB1/$pdb_chain"."_pair_out.txt") or die "$!";    
    while (<INN>) {
        chomp;
        $_ =~ s/^\s+//;
        my @tem = split(/\s+/);
        if (!exists $aa{$tem[0]} and !exists $aa{$tem[3]}) {     
            if (exists $wt_inter_nul{"$tem[4]\t$tem[5]"} and exists $wt_inter_nul{"$tem[1]\t$tem[2]"}) {
                foreach my $j(7..$#tem) {
                    $wt_nul_energy[$j-7] += $tem[$j];   
                }
            }
        }elsif(exists $aa{$tem[0]} and exists $aa{$tem[3]}) {    
            if (exists $wt_inter_pro{$tem[2]} and exists $wt_inter_pro{$tem[5]}) {
                foreach my $j(7..$#tem) {
                    $wt_pro_energy[$j-7] += $tem[$j];
                }
            }
        }else{
                    
            if (exists $wt_inter_nul{"$tem[4]\t$tem[5]"} and exists $wt_inter_pro{$tem[2]}) {
                foreach my $j(7..$#tem) {
                    $wt_inter_energy[$j-7] += $tem[$j];
                }
            }
            if (exists $wt_inter_nul{"$tem[1]\t$tem[2]"} and exists $wt_inter_pro{$tem[5]}) {
                foreach my $j(7..$#tem) {
                    $wt_inter_energy[$j-7] += $tem[$j];
                }
            }

        }
    }
    close INN;
    
    foreach my $k(@wt_inter_energy) {
        print OUT "$k\t";
    }

    foreach my $k(@wt_nul_energy) {
        print OUT "$k\t";
    }
    foreach my $k(@wt_pro_energy) {
        print OUT "$k\t";
    }
    print OUT "\n";
close OUT;

open(OUT,">","$workdir/$mut_name/energy_feature/GB1/GB1_MT_EPI.txt") or die "$!"; 
my %mt_inter_pro;            
my %mt_inter_nul;            
open(INI,"<","$workdir/$mut_name/distance/$pdb_chain"."_pro_nul_dis_5.txt") or die "$!";
while (<INI>) {
        chomp;
        my @tem = split(/\s+/);
        if ($tem[2] <= 5) {
            $mt_inter_pro{$tem[0]} = 0;
            $mt_inter_nul{"$tem[6]\t$tem[5]"} = 0;
        }
    }
    close INI;
    
    my @mt_inter_energy = (0,0,0,0,0);    
    my @mt_pro_energy = (0,0,0,0,0);      
    my @mt_nul_energy = (0,0,0,0,0);      
    
open(INN,"<","$workdir/$mut_name/energy_feature/GB1/$mut_name"."_pair_out.txt") or die "$!";    
while (<INN>) {
        chomp;
        $_ =~ s/^\s+//;
        my @tem = split(/\s+/);
        if (!exists $aa{$tem[0]} and !exists $aa{$tem[3]}) {     
            if (exists $mt_inter_nul{"$tem[4]\t$tem[5]"} and exists $mt_inter_nul{"$tem[1]\t$tem[2]"}) {
                foreach my $j(7..$#tem) {
                    $mt_nul_energy[$j-7] += $tem[$j];   
                }
            }
        }elsif(exists $aa{$tem[0]} and exists $aa{$tem[3]}) {    
            if (exists $mt_inter_pro{$tem[2]} and exists $mt_inter_pro{$tem[5]}) {
                foreach my $j(7..$#tem) {
                    $mt_pro_energy[$j-7] += $tem[$j];
                }
            }
        }else{
                    
            if (exists $mt_inter_nul{"$tem[4]\t$tem[5]"} and exists $mt_inter_pro{$tem[2]}) {
                foreach my $j(7..$#tem) {
                    $mt_inter_energy[$j-7] += $tem[$j];
                }
            }
            if (exists $mt_inter_nul{"$tem[1]\t$tem[2]"} and exists $mt_inter_pro{$tem[5]}) {
                foreach my $j(7..$#tem) {
                    $mt_inter_energy[$j-7] += $tem[$j];
                }
            }

        }
    }
    close INN;
    
    foreach my $k(@mt_inter_energy) {
        print OUT "$k\t";
    }

    foreach my $k(@mt_nul_energy) {
        print OUT "$k\t";
    }
    foreach my $k(@mt_pro_energy) {
        print OUT "$k\t";
    }
    print OUT "\n";
close OUT;

open(IN,"<","$workdir/$mut_name/energy_feature/GB1/GB1_MT_EPI.txt") or die "$!";
chomp(my @mt_line=<IN>);
my @mt=split(/\s+/,$mt_line[0]);
close IN;
    
open(INN,"<","$workdir/$mut_name/energy_feature/GB1/GB1_WT_EPI.txt") or die "$!";
chomp(my @wt_line=<INN>);
my @wt=split(/\s+/,$wt_line[0]);
close INN;
    
open(OUT,">","$workdir/$mut_name/energy_feature/GB1/EPI.txt") or die "$!";
#print OUT
foreach my $kk(0..$#mt) {
my $cha = $mt[$kk] - $wt[$kk];
print OUT $cha,"\t";
}
foreach my $kk(0..$#wt) {
    print OUT $wt[$kk],"\t";
}
close OUT;

################################################Calculate ETOR features########################
open(OUT,">","$workdir/$mut_name/energy_feature/GB1/GB1_WT_ETOR.txt") or die "$!";    
my %wt_distance;
open(INI,"<","$workdir/$mut_name/distance/$pdb_name"."_"."$chain_name"."_diff-distance.txt") or die "$!";      
while (<INI>) {
    chomp;                                                             
    my @tem = split(/\s+/);
        $wt_distance{$mutate_pos}{"$tem[2]\t$tem[0]"} = $tem[3];      
    }
close INI;
    
my @wt_H_bond = (0,0,0,0,0);
my @wt_stacking = (0,0,0,0,0);
my @wt_water_mediated = (0,0,0,0,0);
my @wt_hydrophobic = (0,0,0,0,0);
my @wt_other = (0,0,0,0,0);
        
open(INN,"<","$workdir/$mut_name/energy_feature/GB1/$pdb_chain"."_pair_out.txt") or die "$!";        
while (<INN>) {
    chomp;
        $_ =~ s/^\s+//;
        my @tem = split(/\s+/);                
        if (exists $aa{$tem[3]} and $tem[5] eq $mutate_pos and $tem[2] ne $mutate_pos) {    
            next;    
        }elsif(exists $aa{$tem[0]} and $tem[2] eq $mutate_pos and $tem[5] ne $mutate_pos and exists $aa{$tem[3]}) {  
            
            if ($wt_distance{$tem[2]}{"$tem[4]\t$tem[5]"} <= 3) {
                foreach my $j(7..$#tem) {
                    $wt_H_bond[$j-7] += $tem[$j];                                   
                }
            }elsif($wt_distance{$tem[2]}{"$tem[4]\t$tem[5]"} <= 4) {
                foreach my $j(7..$#tem) {
                    $wt_stacking[$j-7] += $tem[$j];
                }
            }elsif($wt_distance{$tem[2]}{"$tem[4]\t$tem[5]"} <= 5) {
                foreach my $j(7..$#tem) {
                    $wt_water_mediated[$j-7] += $tem[$j];
                }
            }elsif($wt_distance{$tem[2]}{"$tem[4]\t$tem[5]"} <= 6) {
                foreach my $j(7..$#tem) {
                    $wt_hydrophobic[$j-7] += $tem[$j];
                }
            }else{
                foreach my $j(7..$#tem) {
                    $wt_other[$j-7] += $tem[$j];
                }
            }
            
        }elsif($tem[2] eq $mutate_pos and $tem[5] eq $mutate_pos and exists $aa{$tem[0]} and exists $aa{$tem[3]}) {     
            foreach my $j(7..$#tem) {
                print OUT "$tem[$j]\t";                 
            }
        }
        
    }
    close INN;
    
    foreach my $i(@wt_H_bond) {
        print OUT "$i\t";               
    }
    foreach my $i(@wt_stacking) {
        print OUT "$i\t";               
    }
    foreach my $i(@wt_water_mediated) {
        print OUT "$i\t";              
    }
    foreach my $i(@wt_hydrophobic) {
        print OUT "$i\t";               
    }
    foreach my $i(@wt_other) {
        print OUT "$i\t";               
    }
    print OUT "\n";
close IN;
close OUT;       

open(OUT,">","$workdir/$mut_name/energy_feature/GB1/GB1_MT_ETOR.txt") or die "$!";    
my %mt_distance;
open(INI,"<","$workdir/$mut_name/distance/$pdb_name"."_"."$chain_name"."_diff-distance.txt") or die "$!";      
while (<INI>) {
    chomp;                                                             
    my @tem = split(/\s+/);
        $mt_distance{$mutate_pos}{"$tem[2]\t$tem[0]"} = $tem[3];      
    }
close INI;
    
my @mt_H_bond = (0,0,0,0,0);
my @mt_stacking = (0,0,0,0,0);
my @mt_water_mediated = (0,0,0,0,0);
my @mt_hydrophobic = (0,0,0,0,0);
my @mt_other = (0,0,0,0,0);
        
open(INN,"<","$workdir/$mut_name/energy_feature/GB1/$mut_name"."_pair_out.txt") or die "$!";        
while (<INN>) {
    chomp;
        $_ =~ s/^\s+//;
        my @tem = split(/\s+/);                
        if (exists $aa{$tem[3]} and $tem[5] eq $mutate_pos and $tem[2] ne $mutate_pos) {    
            next;    
        }elsif(exists $aa{$tem[0]} and $tem[2] eq $mutate_pos and $tem[5] ne $mutate_pos and exists $aa{$tem[3]}) {  
            
            if ($mt_distance{$tem[2]}{"$tem[4]\t$tem[5]"} <= 3) {
                foreach my $j(7..$#tem) {
                    $mt_H_bond[$j-7] += $tem[$j];                                   
                }
            }elsif($mt_distance{$tem[2]}{"$tem[4]\t$tem[5]"} <= 4) {
                foreach my $j(7..$#tem) {
                    $mt_stacking[$j-7] += $tem[$j];
                }
            }elsif($mt_distance{$tem[2]}{"$tem[4]\t$tem[5]"} <= 5) {
                foreach my $j(7..$#tem) {
                    $mt_water_mediated[$j-7] += $tem[$j];
                }
            }elsif($mt_distance{$tem[2]}{"$tem[4]\t$tem[5]"} <= 6) {
                foreach my $j(7..$#tem) {
                    $mt_hydrophobic[$j-7] += $tem[$j];
                }
            }else{
                foreach my $j(7..$#tem) {
                    $mt_other[$j-7] += $tem[$j];
                }
            }
            
        }elsif($tem[2] eq $mutate_pos and $tem[5] eq $mutate_pos and exists $aa{$tem[0]} and exists $aa{$tem[3]}) {     
            foreach my $j(7..$#tem) {
                print OUT "$tem[$j]\t";                 
            }
        }
        
    }
    close INN;
    
    foreach my $i(@mt_H_bond) {
        print OUT "$i\t";               
    }
    foreach my $i(@mt_stacking) {
        print OUT "$i\t";               
    }
    foreach my $i(@mt_water_mediated) {
        print OUT "$i\t";               
    }
    foreach my $i(@mt_hydrophobic) {
        print OUT "$i\t";               
    }
    foreach my $i(@mt_other) {
        print OUT "$i\t";               
    }
    print OUT "\n";
close IN;
close OUT;

open(IN,"<","$workdir/$mut_name/energy_feature/GB1/GB1_MT_ETOR.txt") or die "$!";
chomp(my @mt_line1=<IN>);
my @mt1=split(/\s+/,$mt_line1[0]);
close IN;
    
open(INN,"<","$workdir/$mut_name/energy_feature/GB1/GB1_WT_ETOR.txt") or die "$!";
chomp(my @wt_line1=<INN>);
my @wt1=split(/\s+/,$wt_line1[0]);
close INN;
    
open(OUTT,">","$workdir/$mut_name/energy_feature/GB1/ETOR.txt") or die "$!";
#print OUT
foreach my $kk(0..$#mt1) {
my $cha = $mt1[$kk] - $wt1[$kk];
print OUTT $cha,"\t";
}
foreach my $kk(0..$#wt1) {
    print OUTT $wt1[$kk],"\t";
}
close OUTT;

#####result.txt####
open(OUT,">","$workdir/$mut_name/result.txt") or die "$!";
print OUT "PDB ID:$pdb_name\nChian ID:$chain_name\nPosition:$mutate_res";
#close OUT;

#########################Random forest regression#############################
system "/home/yaojiang/anaconda3/bin/python $workdir/scripts/MPD-regression.py $workdir";

#########################Random forest classification#############################
system "/home/yaojiang/anaconda3/bin/python $workdir/scripts/MPD-classification.py $workdir";
