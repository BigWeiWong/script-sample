#!/usr/bin/perl 

##############################################################
#  Script     : sample.pl
#  Author     : Qiongzi Qiu
#  Date       : 01/01/2019
#  Last Edited: 01/01/2019, Qiongzi Qiu
#  Description: some description
##############################################################

use strict;
use warnings;

my $num = scalar(@ARGV);
my $usage = "
      Usage: sample.pl file1 file2 file3
        file1:    name of infile1
        file2:    name of infile2
        file3:    name of outfile
";
die "$usage\n" if $num < 1;
die "Errors in arguments\n" if $num != 3;
 
my ($file1, $file2, $file3) = @ARGV;

open(IN1, "$file1") || die "Cannot open the file $file1: $!";
open(IN2, "$file2") || die "Cannot open the file $file2: $!";
open(OUT, ">$file3") || die "Cannot open the file $file3: $!";

my ($tr, $m, %exon, %hash2, @a, $head, $seq1, $seq2, $strand, $length);
my ($k3, $k2, $start3, $start2, $end3, $end2, $oldtr, $num_repeat, $string, $trans_id);
$m = 0;
$oldtr = ""; $line = ""; $seq_len = 0; $seq = ""; $sd = "";

while(<IN1>){
    chomp;
    @a = split;
    $length = $a[4]-$a[3]+1;
    if($a[11] =~ /\"(replace)\"/){
        $tr = $1;
    }
    $strand = $a[6];
    $exon{$m} = {start=>$a[3], end=>$a[4], trans_id=>$tr, strand=>$strand, length=>$length,};
    $m++;
}

{
local $/ = "\n\n";
$c = 0;
while(<IN2>){
    chomp;
    my ($head, $seq1, $seq2) = split(/\n/, $_, 3);
    next unless($head && $seq1 && $seq2);
    $seq2 =~ s/\s+//g;
    @a = split(/\s+/, $head);
    $hash2{$c} = {start=>$a[2], end=>$a[3], seq=>$seq2,};
    $c++;
}
}
close(IN1) || die "Cannot close the file $file1: $!";
close(IN2) || die "Cannot close the file $file2: $!";

for($k3 = 0; $k3 < $m; $k3++){
    $start3 = $exon{$k3}{start};
    $end3 = $exon{$k3}{end};
    $trans_id = $exon{$k3}{trans_id};
    $length = $exon{$k3}{length};
    $strand = $exon{$k3}{strand};
    for($k2 = 0; $k2 < $c; $k2++){
        $start2 = $hash2{$k2}{start}; 
        $end2 = $hash2{$k2}{end};
        $seq2 = $hash2{$k2}{seq};
        next if (($end2 < $start3) || ($start2 > $end3));
        if($end2 <= $end3){
            if($start2 <  $start3){
                $len = $end2 - $start3 + 1;
                $st = $start3 - $start2;
                $string = substr($seq2, $st, $len);
                $seq .= $string;
            }else{
                $seq_len = length($seq);
                $num_repeat = $start2 - $start3 - $seq_len;
                $seq .= "-" x $num_repeat;
                $seq .= $seq2;
            }
        }else{
            if($start2 >= $start3){
                $seq_len = length($seq);
                $num_repeat = $start2 - $start3 - $seq_len;
                $seq .= "-" x $num_repeat;
                $len = $end3 - $start2 + 1;
                $string = substr($seq2, 0, $len);
                $seq .= $string;
            }else{
                $len = $end3 - $start3 + 1;
                $st = $start3 - $start2;
                $string = substr($seq2, $st, $len);
                $seq .= $string;
            }
        }
    }
    $seq_len = length($seq);
    if($seq_len < $length){
        $num_repeat = $length - $seq_len;
        $seq .= "-" x $num_repeat;
    }
    if($trans_id eq $oldtr){
        $line .= $seq;
    }else{
        if($line =~ /\w+/){
            if($sd eq "-"){
                $line = reverse $line;
                $line =~ tr/ACGTacgt/TGCATGCA/;
            }
            print OUT "\>$oldtr\n$line\n";
        }
        $oldtr = $trans_id;
        $line = $seq;
        $sd = $strand;
    }
    $seq = "";
}
if($line =~ /\w+/){
    if($sd eq "-"){
        $line = reverse $line;
        $line =~ tr/ACGTacgt/TGCATGCA/;
    }
    print OUT "\>$oldtr\n$line\n";
}

close(OUT) || die "Cannot close the file $file3: $!";
