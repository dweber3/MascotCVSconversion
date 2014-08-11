use strict;
use warnings;
use Math::BigFloat;

my $S = 'ACDEFGHIK';
my $X = 2;
my @R = @{pepinfo($S, $X)};
#print "Result is @R.\n";
print "Peptide mass is $R[0].\n";
my @distND = @{$R[1]};
print "distND = @distND\n";
print "maxND = $R[2].\n";
print "maxD = $R[3].\n";

sub pepinfo
{
	#pepinfo
	#@_[0] = Scalar or Array Reference subSeq = sequence of amino acids making up the protein
	#@_[1] = Scalar X = exclude N-terminal X residues (default 2)
	#Returns Array Reference R
	#R[0] = Scalar peptideMass
	#R[1] = Array Reference distND = ?
	#R[2] = Scalar maxND  = minimum possible deuterium uptake?
	#R[3] = Scalar maxD = maximum possible deuterium uptake?
	my @R;
	
	
	my $subSeq;
	if (scalar(@_) < 1) {die "No input to parse.\n"};
	if (ref($_[0]) eq '')
	{
		if ($_[0] !~ /^[AC-Z]+$/) {die ('No amino acid sequence passed.\n')}
		else {$subSeq = $_[0];}
	}
	elsif (ref ($_[0]) eq 'ARRAY')
	{
		$subSeq = join ('', @{$_[0]});
	}
	
	my $peptideMass = 0;
	my $C = 0;
	my $N = 0;
	my $O = 0;
	my $S = 0;
	my $Fe = 0;
	my $X = 2;
	my (@distND, $maxND, $maxD, $index);
	
	if (scalar(@_) > 2)
	{$X = $_[1];}
	
	#following values taken from http://en.wikipedia.org/wiki/Proteinogenic_amino_acid
	#heme(C34H31N4O4Fe) from http://www.lfd.uci.edu/~gohlke/molmass/?q=C34H31N4O4Fe
	#acylation(-COCH3)
	
	my @AAshort		= ('A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'Y', '*', ']');
	
	my @AAcarbonNum		= (  3,   3,   4,   5,   9,   2,   6,   6,   6,   6,   5,   4,  12,   5,   5,   6,   3,   4,   3,   5,  11,   9,  34,   2);
	my @AAnitrogenNum	= (  1,   1,   1,   1,   1,   1,   3,   1,   2,   1,   1,   2,   3,   1,   2,   4,   1,   1,   1,   1,   2,   1,   4,   0);
	my @AAoxygenNum		= (  1,   1,   3,   3,   1,   1,   1,   1,   1,   1,   1,   2,   3,   1,   2,   1,   2,   2,   1,   1,   1,   2,   4,   1);
	my @AAsulfurNum		= (  0,   1,   0,   0,   0,   0,   0,   0,   0,   0,   1,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0);
	my @AAironNum		= (  0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   1,   0);
	
	my @AAmonoMass		= (71.03711, 103.00919, 115.02694, 129.04259, 147.06841, 57.02146, 137.05891, 113.08406, 128.09496, 113.08406, 131.04049, 114.04293, 255.15829, 97.05276, 128.05858, 156.10111, 87.03203, 101.04768, 150.95364, 99.06841, 186.07931, 163.06333, 613.1741, 42.0106);
	
	for (my $i = 0; $i < length($subSeq); $i ++)
	{
		print "i = $i.\tAA = " . substr($subSeq, $i, 1) . ".\n\n";
		
		#++$index until $AAshort[$index] eq substr($subSeq, $i); #WARNING! INFINITE LOOP
		for ($index = 0; $index < scalar @AAshort && $AAshort[$index] ne substr($subSeq, $i, 1); $index ++)
		{
			print "index = $index.\tAA_test = $AAshort[$index].\n";
		}
		
		$peptideMass = $peptideMass + $AAmonoMass[$index];
		$C = $C + $AAcarbonNum[$index];
		$N = $N + $AAnitrogenNum[$index];
		$O = $O + $AAoxygenNum[$index];
		$S = $S + $AAsulfurNum[$index];
		$Fe = $Fe + $AAironNum[$index];
	}
	
	$peptideMass = $peptideMass + (1.007825 * 2 + 15.994915); #peptide's mass is the sum of the residue masses plus the mass of water.
	
	###calculate maxND:
	my $obsCThreshold = 1e-3; #set threshold
	
	my $pC13 = 0.0111; #natural richness of C13
	my @CA = (0 .. $C);
	my @distC = @{binopdf(\@CA, $C, $pC13)}; #originally called MATLAB function binopdf(), now calls binopdf() implemented above
	
	
	my $pN15 = 0.00364; #natural richness of N15
	my @NA = (0 .. $N);
	my @distN= @{binopdf(\@NA, $N, $pN15)};

	
	my $pO18 = 0.00205; #natural richness of O18
	my @OA = (0 .. $O);
	my @dist = @{binopdf(\@OA, $O, $pO18)};
	my @distO;
	for (my $i = 0; $i < $O; $i ++)
	{
		@distO[$i*2 - 1] = $dist[$i];
	}
	
	# pS33=0.00762; %natural richness of S33 [ignored here]
	my $pS34=0.04293; #%natural richness of S34
	my @distS;
	if ($S>0)
	{
		my @SA = (0 .. $S);
		@dist = @{binopdf(\@SA,$S,$pS34)};
		for (my $i=0; $i < $S; $i++)
		{
			$index = $i * 2 - 1;
			$distS[$index] = $dist[$i];
		}
    }
	else {@distS = 1;}
	
	my $pFe56=0.91754; #natural richness of Fe56
	my $pFe57=0.02119; #natural richness of Fe57
	my @distFe;
	if ($Fe>0)
	{
		my @FA = (0 .. $Fe);
		@dist = @{binopdf(\@FA,$Fe,$pFe56)}; #//this calc is considering from Fe54 (natural richness 0.05845)
		for (my $i = 1; $i < $Fe; $i++)
		{
			@distFe[$i*2 - 1] = $dist[$i];
		}
		@distFe = @{conv(\@distFe, binopdf(\@FA, $Fe, $pFe57))};
	}
	else
	{
		@distFe=1;
	}
	
	my @finalDist = @{conv(\@distFe, conv(\@distS, conv(\@distO, conv(\@distC, \@distN))))};
	
	$maxND = @finalDist - 1;
	for (my $i = 3; $i < $maxND; $i++)
	{
		if ($finalDist[$i] < $obsCThreshold && $finalDist[$i-1] < $obsCThreshold && $finalDist[$i-2] >= $obsCThreshold)
		{
			$maxND = $i - 3;
			last;
		}
	}
	
	for (my $i = 0; $i < $maxND; $i++)
	{
		$distND[$i] = $finalDist[$i];
	}
	
	#%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	#%%%calculate maxD:
	$maxD = $subSeq - $X;
	for (my $i = $X; $i < length ($subSeq); $i ++)
	{
		if (substr($subSeq, $i) eq 'P')  #exclude Proline
		{
			$maxD--;
		}
	}
	
	@R = ($peptideMass, \@distND, $maxND, $maxD);
	return \@R;
}

sub conv
{
		#Convolution
		#$_[0] = Array Reference K = Kernel
		#$_[1] = Array Reference D = Data
		#Returns Array Reference R
	#print "conv($_[0], $_[1])\n";
	my @K = @{$_[0]};
	#print "Kernel is @K.\n";
	my $K = @K;
	#print "Kernel is $K long.\n";
	my @D = @{$_[1]};
	#print "Data is @D.\n";
	my $D = @D;
	#print "Data is $D long.\n";
	my $R = ($K + $D - 1);
	#print "Result will be $R long.\n";
	my @R = (0 .. ($R - 1));
	for (my $i = 0; $i < $R; $i++)
	{
		#outer loop
		$R[$i] = 0;
		for (my $j = 0; $j <= $i; $j ++)
		{
			#inner loop
			print "\t\t(i, j) are ($i, $j).\n";
			my $Kj = 0;
			if ($j < $K) {$Kj = $K[$j];}
			my $Dij = 0;
			if (($i - $j) < $D) {$Dij = $D[$i - $j];}
			$R[$i] += ($Kj * $Dij);
			print "K[$j] = $Kj. D[$i - $j] = $Dij.\n";
			print "\t R[$i] = $R[$i].\n";
			#print "";
		}
	}
	#print "Result is @R.\n";
	return \@R;
}

sub Abinopdf
{
		#binary probability density function
		#$_[0] = Array Reference X = List of target numbers
		#$_[1] = Scalar N = Number of trials
		#$_[2] = Scalar P = Probability
		#Returns Array Reference R = List of probabilities of target numbers of trials succeeding
	
	#array dereferencing here
	my @X = @{$_[0]};
	my $N = Math::BigFloat->new($_[1]);
	my $P = Math::BigFloat->new($_[2]);
	my $Q = Math::BigFloat->new(1);
	$Q = $Q->bsub($P);
	print "\n\t\tSubroutine Input Array:\n\n\t@X\n";
	my @R = @X;
	
	foreach(@R)
	{
		#do maths on $_
		my $PX = $P->copy()->bpow($_);
		my $NX = $N->bsub($_);
		my $QNX = $Q->bpow($NX);
		$_ = Math::BigFloat->new(($N->bnok($_))->bmul($PX))->bmul($QNX);
	}
	
	print "\n\t\tSubroutine Output Array:\n\n\t@R\n";
	return \@R;
}

sub Sbinopdf
{
		#binomial probability density function
		#$_[0] = Scalar X = Target number
		#$_[1] = Scalar N = Number of trials
		#$_[2] = Scalar P = Individual Probability
		#Returns Scalar R = Probability of target number of trials succeeding
	
	my $X = Math::BigFloat->new($_[0]);
	my $N = Math::BigFloat->new($_[1]);
	my $P = Math::BigFloat->new($_[2]);
	#print "Init: X=$X, N=$N, P=$P\n";
	my $Q = Math::BigFloat->new(1);
	$Q = $Q->bsub($P);
	my $PX = $P->bpow($X);
	#print "P: $P X: $X P^X: $PX\n";
	#print "Q: $Q ";
	my $NX = $N->bsub($X);
	#print "N-X: $NX ";
	my $QNX = $Q->bpow($NX);
	#print "Q^(N-X): $QNX\n\n";
	my $R = Math::BigFloat->new(($N->bnok($X))->bmul($PX))->bmul($QNX);
	return $R;
}