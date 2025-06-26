#!/usr/bin/perl -w
use File::Basename;
use strict;
use threads;
use threads::shared;

unless (@ARGV>=3){
        die "perl $0 <Read1.fastq.gz> <Read2.fastq.gz> <CB_Whilelist.txt> <Thread def:20>\nUsage: Transfer scDR library reads into the fq files listed in the whitelist, then used to run STAR --solo SMART-seq\nNote: 1) only keep the one with correct universal adaptor (pos and seq).\n      2) multi-threads joining and writing the files simultaneously.\nAUTHOR:  Wang Kaile  MDACC llee9050445\@gmal.com\nDATE:  v1:2020/12/28\n";
        }
my ($fq1,$fq2,$wl)=($ARGV[0],$ARGV[1],$ARGV[2]);
my $thr_n;
if (scalar @ARGV == "3"){
$thr_n=20;
}else{
$thr_n=$ARGV[3];
}
my $filename1= basename $fq1;
my $filename2= basename $fq2;
my @file_n = split(/\./,$filename2);
my $folder_pre = $file_n[0];

open (my $FQ1, "gunzip -c $fq1 |") || die "can't open the $fq1 file!";
open (my $FQ2, "gunzip -c $fq2 |") || die "can't open the $fq2 file!";
open(WL, $wl) or die "can't open the $wl file!";
chomp(our @whl=<WL>);
my $folder_name = "${folder_pre}_fq_list_rna_m";
mkdir $folder_name;

#---multi thread
my $total_lines=`zcat $fq1 | wc -l`;
my $thread_num=$thr_n;
my $size = (int($total_lines / ($thread_num * 4)) + 1)*4;

my (%handles1, %handles2);
foreach(1..$thread_num){
    my $outfile1 = $fq1.".$_.tmp";
    my $outfile2 = $fq2.".$_.tmp";
    open $handles1{$_}, ">$outfile1" or die;
    open $handles2{$_}, ">$outfile2" or die;
}

#---split FQ1
my $number1 = 1;
while (<$FQ1>){
    chomp;
    if ($. % $size != 0){
        print {$handles1{$number1}} "$_\n";
    }else{
        print {$handles1{$number1}} "$_\n";
        $number1++;
    }
}
#---split FQ2
my $number2 = 1;
while (<$FQ2>){
    chomp;
    if ($. % $size != 0){
        print {$handles2{$number2}} "$_\n";
    }else{
        print {$handles2{$number2}} "$_\n";
        $number2++;
    }
}

#---start multi-thread
my @fq1_files=glob ("$fq1*.tmp");
my @fq2_files=glob ("$fq2*.tmp");

foreach (0..$thread_num-1){
    #print "Start_thread $_\t";
    my $thr = threads -> create(\&get_func, $fq1_files[$_], $fq2_files[$_]);
}
print "\n";
#---Join results
my %fh2;
for my $label2 (@whl){
my $filename = "${label2}.fq.gz";
open $fh2{$label2}, " | gzip > $folder_name/$filename";
}

while(threads -> list()){
    foreach my $t(threads -> list(threads::joinable)){
	my $tmp = $t -> join();	
	for my $label (@whl){
	if (defined ${tmp} -> {$label}){
	if (defined $fh2{$label}){
	print { $fh2{$label} } ${tmp} -> {$label};
	}
	}
	}
}
}

unlink glob "$fq1*.tmp";
unlink glob "$fq2*.tmp";

#-- for 
sub get_func {
my $fq1s = $_[0];
my $fq2s = $_[1];
open my $fq1sub, $fq1s or die;
open my $fq2sub, $fq2s or die;

my %fh_sub;
foreach my $i (@whl){
$fh_sub{$i} = undef;
}
while(1){
my $a1 = <$fq1sub>;
chomp (my $a2 = <$fq1sub>);
my $a3 = <$fq1sub>;
my $a4 = <$fq1sub>;
my $b1 = <$fq2sub>;
my $b2 = <$fq2sub>;
my $b3 = <$fq2sub>;
my $b4 = <$fq2sub>;
my @c2=split(//,$a2);
my $adp = join("",@c2[18..31]);
my $polyt = join("",@c2[40..47]);
my $cb =join("",@c2[0..7]).join("",@c2[32..39]);

if ($adp eq "GAGGCGTAGTGGCT" && $polyt eq "TTTTTTTT" && "@main::whl" =~ /$cb/){
if (defined $fh_sub{$cb}) {
$fh_sub{$cb} = $fh_sub{$cb}.$b1.$b2.$b3.$b4;
}else{
$fh_sub{$cb} = $b1.$b2.$b3.$b4;
}
}
if(eof($fq1sub))
{
last;
}
}
return \%fh_sub;
}

close $FQ1;
close $FQ2;

