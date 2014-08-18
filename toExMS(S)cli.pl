#		toExMS(S)cli.pl
#	
#	Arranges information in Mascot-exported CVS files into the text files necessary for
#	ExMS. Assumes that scan title includes the retention time in the for "RT:$time".
#	Made to work with Mascot as run by the IBBR as of 2014-07-18. Mascot claims last update 2010-03-30.
#
#	Based on the pepinfo.m function written as part of the ExMS program by ZhongYuan Kan et. al. at the University of Pennsylvania.
#	
###############################################################################
#	Required Export options:
#	(If it isn't listed, this doesn't care about it.)
#
#	Protein Hit Information 				(*)
#
#	Peptide Match Information				(*)
#		Experimental Charge					(*)
#		Calculated Mr (Da)					(*)
#		Mass error (Da)						(*)
#		Start								(*)
#		End									(*)
#		Score								(*)
#		Sequence							(*)
#		Query title							(*)
#
#	Sample line follows:
#	prot_hit_num,prot_acc,pep_query,pep_rank,pep_isbold,pep_isunique,pep_exp_mz,pep_exp_z,pep_calc_mr,pep_delta,pep_start,pep_end,pep_score,pep_res_before,pep_seq,pep_res_after,pep_scan_title
#	1,NLGT,581,1,1,1,520.2807,1,519.2727,0.0007,130,134,31.09,A,SVVCL,L,File:06262014-Fab-MSMS-1.mzXML Scans:801 RT:4.6696min Charge:1+ Fragmentation:cid
###############################################################################
use strict;
use warnings;
use Math::BigFloat;
#open (DEBUG, ">log.txt");

sub binopdf
{
		#binary probability density function
		#$_[0] = Array Reference X = List of target numbers
		#$_[1] = Scalar N = Number of trials
		#$_[2] = Scalar P = Probability
		#Returns Array Reference R = List of probabilities of target numbers of trials succeeding
	
	#array dereferencing here
	my @X = @{$_[0]};
	my $N = new Math::BigFloat $_[1];
	my $P = new Math::BigFloat $_[2];
	my $Q = new Math::BigFloat '1';
	#my $N = Math::BigFloat->new($_[1]);
	#my $P = Math::BigFloat->new($_[2]);
	#my $Q = Math::BigFloat->new(1);
	$Q -= $P;
	#print "\n\t\tSubroutine Input Array:\n\n\t@X\n";
	my @R = @X;
	
	#print "Format: X\tN\tP\tQ\tPX\tNX\tQNX\n\n";
	foreach(@R)
	{
		#do maths on $_
		#printf ("\tRun:\nState = (%6f, %6f, %6f, %6f)\n", $_, $N, $P, $Q);
		#my $PX = $P->bpow($_);
		my $PX = $P ** $_;
		#printf ("PX:\nState = (%6f, %6f, %6f, %6f, %6f)\n", $_, $N, $P, $Q, $PX);
		#my $NX = $N->bsub($_);
		my $NX = $N - $_;
		#printf ("NX:\nState = (%6f, %6f, %6f, %6f, %6f, %6f)\n", $_, $N, $P, $Q, $PX, $NX);
		#my $QNX = $Q->bpow($NX);
		my $QNX = $Q ** $NX;
		#printf ("QNX:\nState = (%6f, %6f, %6f, %6f, %6f, %6f, %6f)\n", $_, $N, $P, $Q, $PX, $NX, $QNX);
		$_ = Math::BigFloat->new(($N->copy()->bnok($_))->bmul($PX))->bmul($QNX);
		#printf ("\n\tResult computed, \$_ = %6f\n\n", $_);
	}
	
	#print "\n\t\tSubroutine Output Array:\n";
	#printf "%6f " x @R . "\n", @R;
	return \@R;
}

sub conv
{
	#Convolution
	#$_[0] = Array Reference K = Kernel
	#$_[1] = Array Reference D = Data
	#Returns Array Reference R
	#print "\t\tconv($_[0], $_[1])\n";
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
			#print "\t\t(i, j) are ($i, $j).\n";
			my $Kj = 0;
			if ($j < $K) {$Kj = $K[$j];}
			my $Dij = 0;
			if (($i - $j) < $D) {$Dij = $D[$i - $j];}
			$R[$i] += ($Kj * $Dij);
			#print "K[$j] = $Kj. D[$i - $j] = $Dij.\n";
			#print "\t R[$i] = $R[$i].\n";
			#print "";
		}
	}
	#print "Result is @R.\n";
	return \@R;
}

sub pepinfo
{
	#pepinfo
	#@_[0] = Scalar or Array Reference subSeq = sequence of amino acids making up the protein
	#@_[1] = Scalar X = exclude N-terminal X residues (default 2)
	#Returns Array Reference R
	#R[0] = Scalar peptideMass
	#R[1] = Array Reference Reference distND = ?
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
	#print "\tpepinfo($subSeq, 2)\n";
	
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
		#print "\ni = $i.\tAA = " . substr($subSeq, $i, 1) . ".\n";
		$index = 0;
		++$index until $AAshort[$index] eq substr($subSeq, $i, 1);
		for ($index = 0; $index < scalar @AAshort && $AAshort[$index] ne substr($subSeq, $i, 1); $index ++)
		{
			#print "index = $index.\tAA_test = $AAshort[$index].\n";
		}
		#print "Matched AA = $AAshort[$index].\n\n";
		
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
	#print "pC13 = $pC13.\n";
	my @CA = (0 .. $C);
	#print "CA:\n@CA\n";
	my @distC = @{binopdf(\@CA, $C, $pC13)}; #originally called MATLAB function binopdf(), now calls binopdf() implemented above
	#print "distC:\n";
	#printf "%6f " x @distC . "\n", @distC;
	
	
	my $pN15 = 0.00364; #natural richness of N15
	my @NA = (0 .. $N);
	my @distN= @{binopdf(\@NA, $N, $pN15)};
	#print "distN:\n@distN\n\n";
	
	my $pO18 = 0.00205; #natural richness of O18
	my @OA = (0 .. $O);
	my @dist = @{binopdf(\@OA, $O, $pO18)};
	my @distO = (0) x ($O * 2);
	#print "distO before:\n@distO\n\n";
	for (my $i = 0; $i < $O; $i ++)
	{
		@distO[($i + 1)*2 - 2] = $dist[$i];
	}
	#print "distO after:\n@distO\n\n";
	
	# pS33=0.00762; %natural richness of S33 [ignored here]
	my $pS34=0.04293; #%natural richness of S34
	my @distS = (0) x ($S * 2);
	if ($S>0)
	{
		my @SA = (0 .. $S);
		@dist = @{binopdf(\@SA,$S,$pS34)};
		for (my $i=0; $i < $S; $i++)
		{
			$index = ($i + 1) * 2 - 2;
			$distS[$index] = $dist[$i];
		}
    }
	else {@distS = 1;}
	#print "distS:\n@distS\n\n";
	
	my $pFe56=0.91754; #natural richness of Fe56
	my $pFe57=0.02119; #natural richness of Fe57
	my @distFe = (0) x ($Fe * 2);
	if ($Fe>0)
	{
		my @FA = (0 .. $Fe);
		@dist = @{binopdf(\@FA,$Fe,$pFe56)}; #//this calc is considering from Fe54 (natural richness 0.05845)
		for (my $i = 1; $i < $Fe; $i++)
		{
			@distFe[($i + 1)*2 - 2] = $dist[$i];
		}
		@distFe = @{conv(\@distFe, binopdf(\@FA, $Fe, $pFe57))};
	}
	else
	{
		@distFe=1;
	}
	#print "distFe:\n@distFe\n\n";
	
	my @interDist = @{conv(\@distC, \@distN)};
	#print "\tconv(distC, distN):\n";
	#printf "%6f " x @interDist . "\n", @interDist;
	
	@interDist = @{conv(\@distO, \@interDist)};
	#print "\tconv(distO, interDist):\n";
	#printf "%6f " x @interDist . "\n", @interDist;
	
	@interDist = @{conv(\@distS, \@interDist)};
	#print "\tconv(distS, interDist):\n";
	#printf "%6f " x @interDist . "\n", @interDist;
	
	@interDist = @{conv(\@distFe, \@interDist)};
	#print "\tconv(distFe, interDist):\n";
	#printf "%6f " x @interDist . "\n", @interDist;
	my @finalDist = @interDist;
	
	#print "finalDist size: " . (scalar @finalDist) . "\n";
	$maxND = (scalar @finalDist) - 1;
	
	#original pepinfo.m loop to identify maxND.
	#for m=3:(maxND+1)
		#if finalDist(m)<obsCThreshold && finalDist(m-1)<obsCThreshold && finalDist(m-2)>=obsCThreshold
			#maxND=m-3; break
		#end
	#end
	for (my $i = 2; $i < $maxND; $i++)
	{
		#print "$obsCThreshold\t> $finalDist[$i]\t> " . $finalDist[$i-1] . "\t< " . $finalDist[$i-2] . "\n";
		if ($finalDist[$i] < $obsCThreshold && $finalDist[$i - 1] < $obsCThreshold && $finalDist[$i - 2] >= $obsCThreshold)
		{
			$maxND = $i - 2;
			last;
		}
	}
	
	for (my $i = 0; $i <= $maxND; $i++)
	{
		$distND[$i] = $finalDist[$i];
	}
	
	#%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	#%%%calculate maxD:
	$maxD = length ($subSeq) - $X;
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

my $infilename = "";
my $outfilename = "";

#get $infilename
print "Please enter the input file:\t";
$infilename = <STDIN>;
chomp $infilename;
#print $DEBUG "Input file is $infilename.\n";
if (-e $infilename && -r $infilename)
{
	open (IFH, '<', $infilename) or die "Could not open input file $infilename.\n";
	#get $outputfilename
	print "Enter desired output filename. If none is entered, file will be `out.txt'.\t";
	$outfilename = <STDIN>;
	chomp $outfilename;
	if ($outfilename eq "") {$outfilename = "out.txt";}
	#print $DEBUG "Output file will be $outfilename.\n";
	open (OFH, '>', $outfilename) or die "Could not open output file $outfilename.\n";
	if (-w $outfilename)
	{
		#do all the things to ifh
		my ($start, $end, $z, $mz, $score, $delta, $seq, $scan, @row, @AoR);
		my $mode = 0;
		my $max = 0;
		#my $i = 0;
		while (<IFH>)
		{
			#do all the things to $_
			if ($_ =~ m/^prot_*+/) #header row
			{
				#assign locations to variables
				$mode = 1;
				#print DEBUG $_;
				@row = split(/,/, $_);
				for (my $i = 0; $i < scalar @row; $i++)
				{
					#ID fields
					if ($row[$i] eq "pep_exp_mz") {$mz = $i;}
					elsif ($row[$i] eq "pep_exp_z") {$z = $i;}
					elsif ($row[$i] eq "pep_delta") {$delta = $i;}
					elsif ($row[$i] eq "pep_score") {$score = $i;}
					elsif ($row[$i] eq "pep_start") {$start = $i;}
					elsif ($row[$i] eq "pep_end") {$end = $i;}
					elsif ($row[$i] eq "pep_seq") {$seq = $i;}
					elsif ($row[$i] eq "pep_scan_title\n") {$scan = $i;}
					else {}
				}
			}
			elsif ($_ =~ m/^\d.+$/ && $mode == 1)	#as long as $_ has data
			{
				#print "Acquiring data row #$i.\n";
				#$i++;
				my $rt;
				@row = split(/,/, $_);
				#call pepinfo
				my @pepinfo = @{pepinfo($row[$seq], 2)};
				if ($row[$scan] =~ /.*RT:(\d+(?:\.\d+)).*/) {$rt = $1;}
				chomp($row[$scan]);
				my $R = $row[$start] . "\t" . $row[$end] . "\t" . $row[$z] . "\t" . $row[$mz] . "\t" . $pepinfo[2] . "\t" . $pepinfo[3] . "\t" . $rt . "\t" . $row[$score] . "\t" . $row[$delta] . "\t";
				foreach (@{$pepinfo[1]})
				{
					$R = $R . $_ . "\t";
				}
				my $columns = $R =~ tr/\t//;
				if ($columns > $max) {$max = $columns;}
				push @AoR, $R;
			}
			else {$mode = 0;} #no more data to read
		}
		#$i = 0;
		foreach (@AoR)
		{
			#print "Outputting row #$i.\n";
			#$i++;
			#pad each row to length
			my $columns = $_ =~ tr/\t//;
			for (my $i = 0; $i < ($max - $columns); $i++)
			{
				$_ = $_ . "0\t";
			}
			chop ($_);
			$_ = $_ . "\n";
			#print @AoR to output here:
			print OFH $_;
		}
	}
	else {die "Output file $outfilename is not writeable.\n";}
}
else
{
	if (!-e $infilename) {die "Input file $infilename does not exist.\n";}
	elsif (!-r $infilename) {die "Input file $infilename is not readable.\n";}
	else {die "This should be unreachable.\n"}
}
#close DEBUG;