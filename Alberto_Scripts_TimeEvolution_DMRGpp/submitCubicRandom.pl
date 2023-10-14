#!/usr/bin/env perl 
use strict;
use warnings;
use utf8;
#use lib '/ccs/home/alnocera/perl5/lib/perl5/ppc64le-linux-thread-multi/';
use Math::Random::MT;
use Math::Trig;

my ($Lx,$Ly,$Lz,$isPeriodic,$J0,$ta,$dt,$m,$seed,$scheduleFile,$precision,$restart) = @ARGV;
my $message = "Lx Ly Lz isPeriodic J0 ta dt m seed scheduleFile.csv precision restart_option\n";
defined($Lx) or die "USAGE: $0 $message";
defined($Ly) or die "USAGE: $0 $message";
defined($Lz) or die "USAGE: $0 $message";
defined($isPeriodic) or die "USAGE: $0 $message";
defined($J0) or die "USAGE: $0 $message";
defined($ta) or die "USAGE: $0 $message";
defined($dt) or die "USAGE: $0 $message";
defined($m) or die "USAGE: $0 $message";
defined($seed) or die "USAGE: $0 $message";
defined($scheduleFile) or die "USAGE: $0 $message";
defined($precision) or die "USAGE: $0 $message";
defined($restart) or die "USAGE: $0 $message";
my $sseedd = $seed+1;
my $gen = Math::Random::MT->new($sseedd);

my $Nsweeps = int($ta*0.5/$dt)+1;
$dt *= pi;

my @matrix;
my $ssize = getSchedule(\@matrix,$scheduleFile);

my $input0 = "in0_Lx$Lx\_Ly$Ly\_J0$J0\_D$m\_seed$seed";

my $input; 

if ($restart eq 'submit' or $restart eq 'test') {
	$input = "in_Lx$Lx\_Ly$Ly\_J0$J0\_ta$ta\_D$m\_seed$seed";
} 
if ($restart eq 'submit_restart' or $restart eq 'test_restart') {
        $input = "inRe_Lx$Lx\_Ly$Ly\_J0$J0\_ta$ta\_D$m\_seed$seed";	
}

my $L = $Lx*$Ly;

my $maxConn = $Ly;

my @matrixConn;

#initialize
for (my $k = 0; $k < $Lz; $k+=1) {
for (my $i = 0; $i < $L; $i+=1) {
	for (my $j = 0; $j < $L; $j+=1) {
		$matrixConn[$i+$L*$j+$L*$L*$k] = 0.0;
	}
}
}
# y direction
for (my $k = 0; $k < $Lz; $k+=1) {
for (my $i = 0; $i < $Lx; $i+=1) {
        for (my $j = 0; $j < $Ly-1; $j+=1) {
                my $x = ($j) + $Ly*($i);
                my $y = ($j+1) + $Ly*($i);
                $matrixConn[$x+$L*$y+$L*$L*$k] = 4.0*JR($J0,$gen,$precision);
		#print STDERR "$x $y $k $matrixConn[$x+$L*$y+$L*$L*$k]\n";
	}
}
}
# x direction
for (my $k = 0; $k < $Lz; $k+=1) {
for (my $i = 0; $i < $Lx-1; $i+=1) {
        for (my $j = 0; $j < $Ly; $j+=1) {
		my $x = ($j) + $Ly*($i);
                my $y = ($j) + $Ly*($i+1);		
                $matrixConn[$x+$L*$y+$L*$L*$k] = 4.0*JR($J0,$gen,$precision);
		#print STDERR "$x $y $k $matrixConn[$x+$L*$y+$L*$L*$k]\n";
	}
}
}
# y periodic
if ($isPeriodic) {
for (my $k = 0; $k < $Lz; $k+=1) {
for (my $i = 0; $i < $Lx; $i+=1) {
                my $x = (0) + $Ly*($i);
                my $y = (0+$Ly-1) + $Ly*($i);
		$matrixConn[$x+$L*$y+$L*$L*$k] = 4.0*JR($J0,$gen,$precision);
		#print STDERR "$x $y $k $matrixConn[$x+$L*$y+$L*$L*$k]\n";		
}	
}
}
# z direction
my @matrixOnsite;
for (my $i = 0; $i < $Lx; $i+=1) {
        for (my $j = 0; $j < $Ly; $j+=1) {
                my $x = ($j) + $Ly*($i);
                my $y = ($j) + $Ly*($i);
                $matrixOnsite[$x+$L*$y] = 4.0*JR($J0,$gen,$precision);
		#print STDERR "$x $y (0-1) $matrixOnsite[$x+$L*$y]\n";
	}
}

createInpo($input0);
createInp($input);

sub JR
{
	my ($J0,$gen,$precision) = @_;
	if ($precision == 1) {
		my $rand = $gen->irand();
        	return $J0*(2.0*($rand%2)-1.0);	
	} else {
                my $rand = $gen->irand();
                my $myrem = $rand % $precision;
		my $mysign = 1;		
		if ($myrem >= ($precision/2)) {
			$mysign =-1;
		}
                my $mymag = ($rand % ($precision/2))+1.0;
                my $myret = $mysign*$mymag*$J0/(1.0*($precision/2));
                return $myret;
	}
}

sub createInpo
{
        my ($input0) = @_;
        open(FOUTIN,"> $input0") or die "$0: Cannot write to $input0 : $!\n";
	my $Lminus2 = $L-2;
	my $Lminus2over2 = ($Lminus2)/2.0;
        my $trunc = 1E-7;
        print FOUTIN<<EOF;
##Ainur 1.0
TotalNumberOfSites=$L;
NumberOfTerms=1;
DegreesOfFreedom=$Lz;
GeometryKind="LongRange";
GeometryOptions="none";
integer GeometryMaxConnections=$maxConn;
matrix Connectors=[
EOF
for (my $i = 0; $i < $L*$Lz; $i+=1) {
        for (my $j = 0; $j < $L*$Lz; $j+=1) {
                if ($j == 0) {print FOUTIN "[";}
                if ($j%$Lz==0 && $i%$Lz==0) {
			my $vv = $matrix[2]*$matrixConn[($i/$Lz)+$L*($j/$Lz)];
                        print FOUTIN "0";			
                } elsif ($j%$Lz==1 && $i%$Lz==1) {
                        my $vv = $matrix[2]*$matrixConn[($i-1)/$Lz+$L*($j-1)/$Lz];			
                        print FOUTIN "0";
                } else {
                        print FOUTIN "0";
                }
                if ($j<$L*$Lz-1) {print FOUTIN ",";}
                if ($j == $L*$Lz-1) {print FOUTIN "]";}
        }
        if ($i<$L*$Lz-1) {print FOUTIN ",\n";}
}
print FOUTIN<<EOF;
];

Model="IsingMultiOrb";
Orbitals=$Lz;

matrix OnSiteLinkSzSz=[
EOF
for (my $k = 0; $k < $Lz-1; $k+=1) {
        for (my $i = 0; $i < $Lx; $i+=1) {
                for (my $j = 0; $j < $Ly; $j+=1) {
                        my $x = ($j)+$Ly*($i);
                        if ($x == 0) {print FOUTIN "[";}
                        my $vv = $matrix[2]*$matrixOnsite[$x+$L*($x)];			
                        print FOUTIN "0";
                        if ($x<$L-1) {print FOUTIN ",";}
                        if ($x == $L-1) {print FOUTIN "]";}
                }
        }
        #if ($k<$Lz-1) {print FOUTIN "],\n[";}  
}
print FOUTIN "];\n";
print FOUTIN<<EOF;
MagneticFieldX=[[
EOF
for (my $k = 0; $k < $Lz; $k+=1) {
for (my $i = 0; $i < $L; $i+=1) {
	my $value = 2.0*$matrix[1];
	print FOUTIN "$value";
	if ($i < $L-1) {print FOUTIN ","}
}
if ($k<$Lz-1) {print FOUTIN "],\n[";}
}
print FOUTIN "]];\n";
print FOUTIN<<EOF;
SolverOptions="twositedmrg,geometryallinsystem,BatchedGemm";
TruncationTolerance=1e-12;
Version="version";
OutputFile="data0_Lz$Lz\_Lx$Lx\_Ly$Ly\_J0$J0\_D$m\_seed$seed";
InfiniteLoopKeptStates=4;
FiniteLoops=[[-$Lminus2, $m, 0],[$Lminus2, $m, 0],[-$Lminus2, $m, 0],[$Lminus2, $m, 0]];

EOF
print STDERR "Written $input0 \n";
close(FOUTIN);
}

sub createInp
{
        my ($input) = @_;
        open(FOUTIN,"> $input") or die "$0: Cannot write to $input: $!\n";
        my $Lminus2 = $L-2;
        my $Lminus2over2 = ($Lminus2)/2.0;
	my $trunc = 1e-10;
	print FOUTIN<<EOF;
##Ainur 1.0
TotalNumberOfSites=$L;
NumberOfTerms=1;
DegreesOfFreedom=$Lz;
GeometryKind="LongRange";
GeometryOptions="none";
integer GeometryMaxConnections=$maxConn;
matrix Connectors=[
EOF
for (my $i = 0; $i < $L*$Lz; $i+=1) {
        for (my $j = 0; $j < $L*$Lz; $j+=1) {
                if ($j == 0) {print FOUTIN "[";}
		if ($j%$Lz==0 && $i%$Lz==0) {
                        my $vv = $matrixConn[($i/$Lz)+$L*($j/$Lz)];			
                	print FOUTIN "$vv";		
		} elsif ($j%$Lz==1 && $i%$Lz==1) {
                        my $vv = $matrixConn[(($i-1)/$Lz)+$L*(($j-1)/$Lz)+$L*$L];			
			print FOUTIN "$vv";                                       
		} else {
			print FOUTIN "0";		
		}
		if ($j<$L*$Lz-1) {print FOUTIN ",";}
                if ($j == $L*$Lz-1) {print FOUTIN "]";}
	}
        if ($i<$L*$Lz-1) {print FOUTIN ",\n";}
}
print FOUTIN<<EOF;
];
Model="IsingMultiOrb";
Orbitals=$Lz;
matrix OnSiteLinkSzSz=[
EOF
for (my $k = 0; $k < $Lz-1; $k+=1) {
	for (my $i = 0; $i < $Lx; $i+=1) {
        	for (my $j = 0; $j < $Ly; $j+=1) {
			my $x = ($j)+$Ly*($i);
                	if ($x == 0) {print FOUTIN "[";}
                        my $vv = $matrixOnsite[($j)+$Ly*($i)+$L*(($j)+$Ly*($i))];			
                	print FOUTIN "$vv";
                	if ($x<$L-1) {print FOUTIN ",";}
                	if ($x == $L-1) {print FOUTIN "]";}
        	}
	}
	#if ($k<$Lz-1) {print FOUTIN "],\n[";}	
}
print FOUTIN "];\n";
print FOUTIN<<EOF;
MagneticFieldX=[[
EOF
for (my $k = 0; $k < $Lz; $k+=1) {
my $value = 2.0;
for (my $i = 0; $i < $L; $i+=1) {
print FOUTIN "$value";
if ($i < $L-1) {print FOUTIN ","}
}
if ($k<$Lz-1) {print FOUTIN "],\n[";}
}
print FOUTIN "]];\n";
print FOUTIN<<EOF;

string PrintHamiltonianAverage="s==c";
SolverOptions="twositedmrg,BatchedGemm,calcAndPrintEntropies,restart,TimeStepTargeting,normalizeTimeVectors,recoveryenableread,noclobber";
RestartFilename="data0_Lz$Lz\_Lx$Lx\_Ly$Ly\_J0$J0\_D$m\_seed$seed";
string RecoverySave="\@M=2,%l%%1";
Version="version";
OutputFile="data1_Lz$Lz\_Lx$Lx\_Ly$Ly\_J0$J0\_D$m\_seed$seed";
InfiniteLoopKeptStates=4;
FiniteLoops=[
EOF
for (my $i = 0; $i < $Nsweeps-1; $i++) {
print FOUTIN "[-$Lminus2, $m, 2],[ $Lminus2, $m, 2],"
}
print FOUTIN<<EOF;
[-$Lminus2, $m, 2],[ $Lminus2, $m, 3]
];
TruncationTolerance=$trunc;
real ta=$ta;
TSPTau=$dt;
TSPTimeSteps=4;
TSPAdvanceEach=$Lminus2;
TSPAlgorithm="Krylov";
TSPSites=[$Lminus2over2];
TSPLoops=[0];
TSPProductOrSum="product";
GsWeight=0.1;
real TridiagEps=1e-8;
int TridiagSteps=200;
string TSPOp0:TSPOperator="expression";
string TSPOp0:OperatorExpression="identity";
Threads=1;

EOF
print FOUTIN "TimeSchedule=[";
for (my $ss = 0; $ss<$ssize-1; $ss++) {
print FOUTIN "[$matrix[3*$ss],$matrix[1+3*$ss],$matrix[2+3*$ss]],\n";
}
print FOUTIN "[$matrix[3*($ssize-1)],$matrix[1+3*($ssize-1)],$matrix[2+3*($ssize-1)]]];\n\n";

print STDERR "Written $input\n";
close(FOUTIN);
}


sub getSchedule
{
        my ($matrix,$scheduleF) = @_;
	open(FIN,$scheduleF) or die "$0: Cannot open $scheduleF: $!\n";
	my $counter = 0;
	while (<FIN>) {
	$counter++;
	}	
	close(FIN);
	my $ssize1 = $counter;
	open(FIN1,$scheduleF) or die "$0: Cannot open $scheduleF: $!\n";
	my $s = 0;
	while (<FIN1>) {
	$_ =~ s/\(//g;
	$_ =~ s/,/ /g;
	$_ =~ s/\)//g;
	chomp;
	my @temp = split;
	@$matrix[3*$s] = $temp[0]; 
        @$matrix[1+3*$s]= -$temp[1]; 
        @$matrix[2+3*$s]= $temp[2]; 	
	$s++;
	}
	close(FIN1);
	return $ssize1;
}
