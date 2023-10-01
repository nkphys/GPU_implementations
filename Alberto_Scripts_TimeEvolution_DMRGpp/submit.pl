#!/usr/bin/env perl 
use strict;
use warnings;
use utf8;

my ($Lx,$Ly,$Lz,$J0,$dt,$isPeriodic,$sche,$precision,$restart,$times) = @ARGV;
my $message = "Lx Ly Lz J0 dt isperiodic? schedule.csv precision [test,submit,submit_restart]\n";
defined($Lx) or die "USAGE: $0 $message";
defined($Ly) or die "USAGE: $0 $message";
defined($Lz) or die "USAGE: $0 $message";
defined($J0) or die "USAGE: $0 $message";
defined($dt) or die "USAGE: $0 $message";
defined($isPeriodic) or die "USAGE: $0 $message";
defined($sche) or die "USAGE: $0 $message";
defined($precision) or die "USAGE: $0 $message";
defined($restart) or die "USAGE: $0 $message";

#my @ta = (1.0,2.0,4.0,6.0,8.0,10.0,12.0,14.0,16.0,18.0,20.0,24.0,28.0,32.0,36.0,40.0,44.0,48.0,50.0);
my @m = (200,250,350,500,750);
#my @ta = (6,10,20,5,8,12,14,16,18,22);
my @ta = (7);#,6,8,10,9,12,14,16,18,20,23,26,30,34);
#my @m = (48);
#
#
my $recoveryfile; 
my $exe = "dmrg";
my $exeO = "observe";
# USAGE: submitSquareRandom.pl Lx Ly isPeriodic J0 ta dt m seed scheduleFile.csv
# create inputs
for (my $k =0; $k<scalar(@ta);$k++) {
	for (my $j =0; $j<scalar(@m);$j++) {

		my $input0 = "in0_Lx$Lx\_Ly$Ly\_J0$J0\_D$m[$j]\_seed\${LSB_JOBINDEX}";
                my $input = "in_Lx$Lx\_Ly$Ly\_J0$J0\_ta$ta[$k]\_D$m[$j]\_seed\${LSB_JOBINDEX}";
		my $inputRe = "inRe_Lx$Lx\_Ly$Ly\_J0$J0\_ta$ta[$k]\_D$m[$j]\_seed\${LSB_JOBINDEX}";
		my $out0 = "out0_Lx$Lx\_Ly$Ly\_J0$J0\_D$m[$j]\_seed\${LSB_JOBINDEX}";
                my $out = "out_Lx$Lx\_Ly$Ly\_J0$J0\_ta$ta[$k]\_D$m[$j]\_seed\${LSB_JOBINDEX}";
                my $out1 = "out1_Lx$Lx\_Ly$Ly\_J0$J0\_ta$ta[$k]\_D$m[$j]\_seed\${LSB_JOBINDEX}";

		if ($restart eq 'submit' or $restart eq 'test') {
                	for (my $i =1; $i<2;$i+=10) {
				system("perl submitCubicRandom.pl $Lx $Ly $Lz $isPeriodic $J0 $ta[$k] $dt $m[$j] $i $sche $precision $restart");
				$recoveryfile = "Recovery2data1_Lz$Lz\_Lx$Lx\_Ly$Ly\_J0$J0\_D$m[$j]\_seed$i.hd5";
                	}
                                my $Batch = "Batch\_ta$ta[$k]\_D$m[$j].pbs";
                                createBatch($Batch,0,$input0,$input,$inputRe,$out0,$out,$out1,$sche,$ta[$k],$J0,$m[$j],$exe,$restart);
			
			if ($restart eq 'submit') {
				system("bsub $Batch");
			}
		} 
                if ($restart eq 'submit_restart' or $restart eq 'test_restart') {

                        my $kk = 0;
                        my $fsub;
                        my $id = 0;			
                        for (my $i =1; $i<2;$i+=10) {
                                system("perl submitCubicRandom.pl $Lx $Ly $Lz $isPeriodic $J0 $ta[$k] $dt $m[$j] $i $sche $precision $restart");
                                $recoveryfile = "Recovery2data1_Lz$Lz\_Lx$Lx\_Ly$Ly\_J0$J0\_D$m[$j]\_seed$i.hd5";
                        }			
			defined($times) or die "USAGE: $0 PROVIDE \"times\" you want to repeat the job!\n";
			my $Batch = "Batch$kk\_ta$ta[$k]\_D$m[$j].pbs";
                        createBatch($Batch,$kk,$input0,$input,$inputRe,$out0,$out,$out1,$sche,$ta[$k],$J0,$m[$j],$exe,$restart);

			if ($restart eq 'submit_restart') {
				$fsub = `bsub $Batch`;
				$id = getid($fsub);
			}
			for ($kk = 1; $kk<=$times; $kk++) {
                                $Batch = "Batch$kk\_ta$ta[$k]\_D$m[$j].pbs";
                                createBatch($Batch,$kk,$input0,$input,$inputRe,$out0,$out,$out1,$sche,$ta[$k],$J0,$m[$j],$exe,$restart);				
				if ($restart eq 'submit_restart') {
					$fsub = `bsub -w ended\\($id\\) $Batch`;
					$id = getid($fsub);
				}
			}
                }
	}
}

sub getid
{
	my ($string) = @_;
	$string =~ s/<//g;
        $string =~ s/>//g;
	my @temp = split(' ', $string);
	print STDERR $string;
	return $temp[1];
}
sub createBatch
{
        my ($Batch,$kk,$input0,$input,$inputRe,$out0,$out,$out1,$sche,$ta,$J0,$m,$exe,$restart) = @_;
        open(FOUTINB,"> $Batch") or die "$0: Cannot write to $Batch : $!\n";
        print FOUTINB <<EOF;
#!/bin/bash
#BSUB -J Lx$Lx\_Ly$Ly\_ta$ta\_D$m.[1-2:10]
#BSUB -P CPH149
#BSUB -W 24:00
#BSUB -nnodes 1
#BSUB -q killable
#BSUB -alloc_flags "gpumps smt1"
#BSUB -o Lx$Lx\_Ly$Ly\_ta$ta\_D$m.%J
#BSUB -e Lx$Lx\_Ly$Ly\_ta$ta\_D$m.%J

module load gcc/9.3.0
module load hdf5/1.10.7
module load cuda/11.0.3
module load openblas
module load magma/2.6.2
module load boost/1.76.0-global_symbols
EOF
if ($kk==0) {
print FOUTINB <<EOF;	
mkdir -p \$MEMBERWORK/cph149/glass_3d_6x6x2_D$m
EOF
}
if ($kk>=0) {
print FOUTINB <<EOF;
cd \$MEMBERWORK/cph149/glass_3d_6x6x2_D$m
EOF
}
if ($kk==0) {
print FOUTINB <<EOF;
cp \$PROJWORK/cph149/glass_3d/6x6x2_openblas_restart/{$sche,$input0,$inputRe,$exe,$exeO} .
date
jsrun -n 1 -r 1 -a 1 -c 7 -bpacked:7 -g 1 ./$exe -f $input0 -p 12 -l $out0
mkdir /ccs/proj/cph149/data_glass_3d/6x6x2/
cp $out0 /ccs/proj/cph149/data_glass_3d/6x6x2/
date
EOF
} 
if ($kk>0){
print FOUTINB <<EOF;
rm $recoveryfile
#cp /ccs/proj/cph149/data_glass_3d/6x6x2/$recoveryfile .
#cp /ccs/proj/cph149/data_glass_3d/6x6x2/$out .
EOF
}
if ($kk>=0) {
print FOUTINB <<EOF;
jsrun -n 1 -r 1 -a 1 -c 7 -bpacked:7 -g 1 ./$exe -f $inputRe -p 12 -l $out "<P0|sx|P0>"
#cp data1*hd5 /ccs/proj/cph149/data_glass_3d/6x6x2/
cp $out /ccs/proj/cph149/data_glass_3d/6x6x2/
EOF
}
if ($kk==$times) {
print FOUTINB <<EOF;
#cp $recoveryfile /ccs/proj/cph149/data_glass_3d/6x6x2/
date
jsrun -n 1 -r 1 -a 1 -c 7 -bpacked:7 -g 1 ./observe -f $inputRe "<P0|sz;sz|P0>" &> $out1
date
cp $out1 /ccs/proj/cph149/data_glass_3d/6x6x2/
rm -r \$MEMBERWORK/cph149/glass_3d_6x6x2_D$m
EOF
}
print STDERR "Written $Batch\n";
}

