use warnings;
use strict;
BEGIN{
     unshift (@INC,'/root/perl5/lib/perl5/5.16.3/x86_64-linux-thread-multi');
     unshift (@INC,'/root/perl5/lib/perl5/5.16.3');
     unshift (@INC,'/root/perl5/lib/perl5/x86_64-linux-thread-multi');
     unshift (@INC,'/root/perl5/lib/perl5');
}
use PDL;
#no warnings 'experimental::smartmatch';
my @aliphatic=(0,19,10,9,12,4);
my @aromatic=(13,17,18,8);
my @polar=(15,16,2,5);
my @positive=(11,1);
my @negative =(3,6);
my @special=(7,14);
my $work="$ARGV[0]/output/bio_feature/JSD";
my @blosums=();

my %hash_blosum=();
my $line_n=0;
open(FA, "<","$ARGV[0]/scripts/BLOSUM62") or die "$!";
while (<FA>)
{
     chomp;
     s/^\s+//g;
     my @item=split/\s+/;
     last if $item[0]=~/^B/;
     #print scalar @item,"\n";
     if (scalar @item>=26)
     {
        my @ref=@item[1..20];
        push(@blosums,[@ref]);
        for(my $i=0;$i<20;$i++)
        {
          $hash_blosum{$line_n}{$i}=$ref[$i];
        }
        $line_n++;
     }
}
close FA;
#print $hash_blosum{0}{0},"\tcheck\n";
my $blosum=pdl(@blosums[0..19]);
#print $blosum,"\n";

opendir(DIR,"$work/pssm");
my @dird = readdir(DIR);
closedir DIR;
my @chains;
foreach my $d(@dird) {
     if ($d =~ /\.pssm/) {
          $d =~ s/\.pssm//;
          push @chains,$d;
     }
}

foreach my $chain (@chains)
{
     if (-e "$work/cscores/$chain.txt") {
          next;
     }else{
     
     #print $chain,"\n";
     $chain=~s/\s+//g;
     my $pssm_file="$work/pssm/$chain.pssm";
     my ($seq,$aa,$se,$re,$dc,$vnes,$spcp,$sl,$win_se,$win_re,$win_dc,$win_vne,$win_spc,$win_sl)=&JSD($pssm_file);
     my @seq_num=@{$seq};
     my @aa_key=@{$aa};
     my @se=@{$se};
     my @re=@{$re};
     my @dc=@{$dc};
     my @vnes=@{$vnes};
     my @spcs=@{$spcp};
     my @sls=@{$sl};
     my @wins_se=@{$win_se};
     my @wins_re=@{$win_re};
     my @wins_dc=@{$win_dc};
     my @wins_vne=@{$win_vne};
     my @wins_spc=@{$win_spc};
     my @wins_sl=@{$win_sl};
     
     my $sen=&zscore_logistic(@se);
     my $ren=&zscore_logistic(@re);
     my $dcn=&zscore_logistic(@dc);
     my $vnesn=&zscore_logistic(@vnes);
     my $spcsn=&zscore_logistic(@spcs);
     my $slsn=&zscore_logistic(@sls);
     my $winsn_se=&zscore_logistic(@wins_se);
     my $winsn_re=&zscore_logistic(@wins_re);
     my $winsn_dc=&zscore_logistic(@wins_dc);
     my $winsn_vne=&zscore_logistic(@wins_vne);
     my $winsn_spc=&zscore_logistic(@wins_spc);
     my $winsn_sl=&zscore_logistic(@wins_sl);
     
     my @sen=@{$sen};
     my @ren=@{$ren};
     my @dcn=@{$dcn};
     my @vnesn=@{$vnesn};
     my @spcsn=@{$spcsn};
     my @slsn=@{$slsn};
     my @winsn_se=@{$winsn_se};
     my @winsn_re=@{$winsn_re};
     my @winsn_dc=@{$winsn_dc};
     my @winsn_vne=@{$winsn_vne};
     my @winsn_spc=@{$winsn_spc};
     my @winsn_sl=@{$winsn_sl};
     
     my $out_cs="$work/cscores/$chain.txt";
     open(OUT, ">$out_cs") or die "$!";
     for(my $i=0;$i<@sen;$i++)
     {
         print OUT $seq_num[$i],"\t",$aa_key[$i],"\t",$sen[$i],"\t",$ren[$i],"\t",$dcn[$i],"\t$vnesn[$i]\t$spcsn[$i]\t$slsn[$i]\t$winsn_se[$i]\t$winsn_re[$i]\t$winsn_dc[$i]\t$winsn_vne[$i]\t$winsn_spc[$i]\t$winsn_sl[$i]\t\n";
     }
     }
}

sub JSD
{
    my($dir_in)=@_;
    my $chain_len=0;
    my @chain=();
    my @seq_num=();
    my @aa_key=();
    my @SE=();
    my @RE=();
    my @DC=();
    my @VNE=();
    my @SPC=();
    my @SL=();
    my @residue_back=(0.078,0.051,0.041,0.052,0.024,0.034,0.059,0.083,0.025,0.062,0.092,0.056,0.024,0.044,0.043,0.059,0.055,0.014,0.034,0.072);
    open(IN,$dir_in)||die "Can not open file: $dir_in in the read_PSSM package\n";
    my @pssm_raw=();
    while(my $raw=<IN>)
    {
          my @raw_wop=();
          my ($se,$re,$dc);
          $raw=~s/\n|\r//g;
          my @raw_inf=split(/ +/,$raw);
          my $num=scalar(@raw_inf);
          my $t=0;
          my $sum=0;
          my ($alip_num,$arom_num,$pos_num,$neg_num,$polar_num,$spec_num);
          if($num>42)
          {
               
               $chain[$chain_len]=$raw_inf[1]."\t".$raw_inf[2];
               $chain_len++;
               for(my $i=0; $i<20;$i++)
               {
                    $raw_wop[$t]=$raw_inf[$i+23];
                    $t++;
               }
               my @vne_matrix=();
               for(my $i=0;$i<20;$i++)
               {
                    $sum+=$raw_wop[$i];
                    my @matrix_initial=(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0);
                    $matrix_initial[$i]=$raw_wop[$i];
                    #print $matrix_initial[$i],"\n";
                    push(@vne_matrix,[@matrix_initial]);
                    if ($i ~~ @aliphatic)
                    {
                         $alip_num+=$raw_wop[$i];
                    }elsif($i ~~ @aromatic)
                    {
                         $arom_num+=$raw_wop[$i];
                    }elsif($i ~~ @polar)
                    {
                         $polar_num+=$raw_wop[$i];
                    }elsif($i~~ @positive)
                    {
                         $pos_num+=$raw_wop[$i];
                    }elsif($i~~ @negative)
                    {
                         $neg_num+=$raw_wop[$i];
                    }elsif($i~~ @special)
                    {
                         $spec_num+=$raw_wop[$i];
                    }else{
                         die;
                    }
                    
                    
               }
               if ($alip_num ==0) {
                    $alip_num =1;
               }
               if ($arom_num ==0) {
                    $arom_num =1;
               }
               if ($polar_num ==0) {
                    $polar_num =1;
               }
               if ($pos_num ==0) {
                    $pos_num =1;
               }
               if ($neg_num ==0) {
                    $neg_num =1;
               }
               if ($spec_num ==0) {
                    $spec_num =1;
               }
               my $sl=-($alip_num* log($alip_num)+$arom_num* log($arom_num)+$polar_num* log($polar_num)+$pos_num* log($pos_num)+$neg_num* log($neg_num)+$spec_num* log($spec_num));
               #vne_matrix
               my $vne=pdl(@vne_matrix[0..19]);
               my $vne_product=($vne x $blosum);
               
               my @plog=();
               for(my $i=0; $i<20;$i++)
               {
                    my @matrix_initial=(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0);
                    my $temp=($vne_product->slice("$i","$i"));
                    if($temp=~/([0-9]+)/g)
                    {
                         #print $1,"\n";
                         if ($1 == 0)
                         {
                              $matrix_initial[$i]=0;
                         }else{
                             $matrix_initial[$i]=&log20($1); #?
                         }
                         push(@plog,[@matrix_initial]);
                    }else{
                         die;
                    }
               }
               
               my @p=();
               for(my $i=0; $i<20;$i++)
               {
                   my @matrix_initial=(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0);
                   my $temp=$vne_product->slice("$i","$i");
                    if($temp=~/([0-9]+)/g)
                    {
                         #print $1,"\n";
                         $matrix_initial[$i]=$1;
                         push(@p,[@matrix_initial]);
                    }else{
                         die;
                    }
               }
               my $plog=pdl(@plog[0..19]);
               my $p=pdl(@p[0..19]);
               my $vne_results=($p x $plog);
               
               my $vne_tr=0;
               for(my $i=0; $i<20;$i++)
               {
                    my $matrix_tr=$vne_results->slice("$i","$i");
                    if($matrix_tr=~/([0-9]+)/g)
                    {
                         $vne_tr+=(-$1);
                    }else{
                         die;
                    }
               }
               
               for(my $i=0;$i<20;$i++)
               {
                    if($sum==0 or $raw_wop[$i]==0)
                    {
                         $raw_wop[$i]=0;
                         $se+=0;
                         $re+=0;
                         $dc+=0.5*$residue_back[$i]* log($residue_back[$i]/(0.5*$raw_wop[$i]+0.5*$residue_back[$i]));
                    }else
                    {
                         $raw_wop[$i]=$raw_wop[$i]/$sum;
                         $se+=-$raw_wop[$i]* log($raw_wop[$i]);
                         $re+=$raw_wop[$i]* log($raw_wop[$i]/$residue_back[$i]);
                         $dc+=0.5*$raw_wop[$i]* log($raw_wop[$i]/(0.5*$raw_wop[$i]+0.5*$residue_back[$i]))+0.5*$residue_back[$i]* log($residue_back[$i]/(0.5*$raw_wop[$i]+0.5*$residue_back[$i]));
                    }			
               }
               my $spc_line=&SPC(@raw_wop);
               push(@seq_num,$raw_inf[1]);
               push(@aa_key,$raw_inf[2]);
               push(@SE,$se);
               push(@RE,$re);
               push(@DC,$dc);
               push(@VNE,$vne_tr);
               push(@SPC,$spc_line);
               push(@SL,$sl);
        }
        
    }
    close(IN);
    my $win_se=&WIN(@SE);
    my $win_re=&WIN(@RE);
    my $win_dc=&WIN(@DC);
    my $win_vne=&WIN(@VNE);
    my $win_spc=&WIN(@SPC);
    my $win_sl=&WIN(@SL);
    return (\@seq_num,\@aa_key,\@SE,\@RE,\@DC,\@VNE,\@SPC,\@SL,$win_se,$win_re,$win_dc,$win_vne,$win_spc,$win_sl);
}

sub zscore_logistic
{
    my @tempsub=@_;
    my $total=0;
    foreach my $item (@tempsub)
    {
        $total+=$item;
    }
    my $avgsub=$total/scalar @tempsub;
    
    my $sdsub=0;
    foreach my $item (@tempsub){
        $sdsub+=($item-$avgsub)**2;
    }
    my $stdsub=sqrt($sdsub/($#tempsub+1));
   # print $avgsub,"\t$stdsub\n";
    my @norm=();
    foreach my $item (@tempsub)
    {
        my $norm=0;
        if ($stdsub ==0)
        {
            $norm=0;
        }else{
            $norm=($item-$avgsub)/$stdsub;
            $norm=sprintf("%.4f",$norm);
        }
        push(@norm,$norm);
    }
    return \@norm;
}
sub log20 {
  my $n = shift;
  return log($n)/log(20);
}

sub SPC
{
     my @pssm_line=@_;
     my $spc_fore=0;
     my $spc_back=0;
     for(my $k=0;$k<20;$k++)
     {
          for(my $q=$k+1;$q<20;$q++)
          {
               if (exists $hash_blosum{$k}{$q})
               {
                    $spc_fore+=$pssm_line[$k]*$pssm_line[$q];
                    $spc_back+=$pssm_line[$k]*$pssm_line[$q]*$hash_blosum{$k}{$q};
                    #print $spc_fore,"\t$spc_back\n";
               }else{
                    #print $k,"\t$q\n";
                    die;
               }
               
          }
     }
     my $spc=1-$spc_fore*$spc_back;
     return $spc;
     
}

sub WIN
{
    #window_3
    my @scores=@_;
    my @WINS=();
    #for(my $i=0;$i<@scores;$i++)
    #{
    #    my $dd=$scores[$i];
    #    if($i==0)
    #    {
    #       my $si_win=$dd+$scores[$i+1];
    #       my $win=0.5*$dd+0.5*$si_win/3;
    #       push(@WINS,$win);
    #    }elsif($i==@scores-1)
    #    {
    #       my $si_win=$scores[$i-1]+$dd;
    #       my $win=0.5*$dd+0.5*$si_win/3;
    #       push(@WINS,$win);
    #    }else{
    #       my $si_win=$scores[$i-1]+$dd+$scores[$i+1];
    #       my $win=0.5*$dd+0.5*$si_win/3;
    #       push(@WINS,$win);
    #    }
    #}
    #return \@WINS;
    
    #================window eq 6======================#
    for(my $i=0;$i <= $#scores;$i++) {
     if ($i eq 0) {
          my $si_win=$scores[$i] + $scores[$i+1] + $scores[$i+2] + $scores[$i+3];
          my $win = 0.5*$scores[$i] + 0.5*($si_win/7);
          push(@WINS,$win);
     }elsif($i eq 1) {
          my $si_win=$scores[$i-1] + $scores[$i] + $scores[$i+1] + $scores[$i+2] + $scores[$i+3];
          my $win = 0.5*$scores[$i] + 0.5*($si_win/7);
          push(@WINS,$win);
     }elsif($i eq 2) {
          my $si_win=$scores[$i-2] + $scores[$i-1] + $scores[$i] + $scores[$i+1] + $scores[$i+2] + $scores[$i+3];
          my $win = 0.5*$scores[$i] + 0.5*($si_win/7);
          push(@WINS,$win);          
     }elsif($i eq $#scores) {
          my $si_win=$scores[$i-3] + $scores[$i-2] + $scores[$i-1] + $scores[$i];
          my $win = 0.5*$scores[$i] + 0.5*($si_win/7);
          push(@WINS,$win);            
     }elsif($i eq ($#scores - 1) ) {
          my $si_win=$scores[$i-3] + $scores[$i-2] + $scores[$i-1] + $scores[$i] + $scores[$i+1];
          my $win = 0.5*$scores[$i] + 0.5*($si_win/7);
          push(@WINS,$win);          
     }elsif($i eq ($#scores - 2)) {
          my $si_win=$scores[$i-3] + $scores[$i-2] + $scores[$i-1] + $scores[$i] + $scores[$i+1] + $scores[$i+2];
          my $win = 0.5*$scores[$i] + 0.5*($si_win/7);
          push(@WINS,$win);           
     }else{
          my $si_win=$scores[$i-3] + $scores[$i-2] + $scores[$i-1] + $scores[$i] + $scores[$i+1] + $scores[$i+2] + $scores[$i+3];
          my $win = 0.5*$scores[$i] + 0.5*($si_win/7);
          push(@WINS,$win);          
     }
    }
    return \@WINS;
}

