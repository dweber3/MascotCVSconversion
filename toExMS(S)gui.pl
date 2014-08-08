#		toExMS(S)gui.pl
#	
#	Arranges information in Mascot-exported CVS files to that necessary for
#	ExMS. Remember to export to .xlsx afterwards! Made to work with Mascot as run by the IBBR as of 2014-07-18.
#	Mascot claims last update 2010-03-30.
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
use Tk;
use strict;
use warnings;
use Math::BigFloat;
open (DEBUG, ">log.txt");

sub binopdf
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

sub pepinfo
{
	if (scalar(@_) < 1) {die "No input to parse.\n"};
	if (@_[0] !~ /^[AC-Z]+$/) {die "No amino acid sequence passed.\n"}
	
	my $subSeq = $_[0];
	my $peptideMass = 0, $distND, $maxND, $maxD, $C = 0, $N = 0, $O = 0, $S = 0, $Fe = 0, $X = 2, $index, @distO;
	
	#following values taken from http://en.wikipedia.org/wiki/Proteinogenic_amino_acid
	#heme(C34H31N4O4Fe) from http://www.lfd.uci.edu/~gohlke/molmass/?q=C34H31N4O4Fe
	#acylation(-COCH3)
	
	@AAshort		= ('A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'Y', '*', ']');
	
	@AAcarbonNum	= (  3,   3,   4,   5,   9,   2,   6,   6,   6,   6,   5,   4,  12,   5,   5,   6,   3,   4,   3,   5,  11,   9,  34,   2);
	@AAnitrogenNum	= (  1,   1,   1,   1,   1,   1,   3,   1,   2,   1,   1,   2,   3,   1,   2,   4,   1,   1,   1,   1,   2,   1,   4,   0);
	@AAoxygenNum	= (  1,   1,   3,   3,   1,   1,   1,   1,   1,   1,   1,   2,   3,   1,   2,   1,   2,   2,   1,   1,   1,   2,   4,   1);
	@AAsulfurNum	= (  0,   1,   0,   0,   0,   0,   0,   0,   0,   0,   1,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0);
	@AAironNum		= (  0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   1,   0);
	
	@AAmonoMass		= (71.03711, 103.00919, 115.02694, 129.04259, 147.06841, 57.02146, 137.05891, 113.08406, 128.09496, 113.08406, 131.04049, 114.04293, 255.15829, 97.05276, 128.05858, 156.10111, 87.03203, 101.04768, 150.95364, 99.06841, 186.07931, 163.06333, 613.1741, 42.0106);
	
	for (my $i = 0; $i < $#subSeq; $i ++)
	{
		++$index until $AAshort[$index] eq substr($subSeq, $i);
		$peptideMass = $peptideMass + $AAmonoMass[$index];
		$C = $C + $AAcarbonNum[$index];
		$N = $N + $AAnitrogenNum[$index];
		$O = $O + $AAoxygenNum[$index];
		$S = $S + $AAsulferNum[$index];
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
		@distO[$i *2 - 1] = $dist[$i];
	}
	
	# pS33=0.00762; %natural richness of S33 [ignored here]
	my $pS34=0.04293; %natural richness of S34
	my @distS;
	if ($S>0)
	{
		my @SA = (0 .. $S);
		@dist = @{binopdf(\@SA,$S,$pS34)};
		for (my $i=0; $i < $S; $i++)
		{
			@distS[$i*2-1]=@dist($i);
		}
    }
	else {@distS=1;}
	
	my $pFe56=0.91754; #natural richness of Fe56
	my $pFe57=0.02119; #natural richness of Fe57
	my @distFe;
	if (Fe>0)
	{
		my @FA = (0..$Fe)
		@dist = @{binopdf(\@FA,$Fe,$pFe56)}; #//this calc is considering from Fe54 (natural richness 0.05845)
		for (my $i = 1; $i < $Fe; $i++)
		{
			@distFe[i*2-1]=@dist[$i];
		}
		@distFe = @{conv(\@distFe, binopdf(\@FA, $Fe, $pFe57))};
	}
	else
	{
		@distFe=1;
	}
	
	my @finalDist = @{conv(\@distFe, conv(\@distS, conv(\@distO, conv(\@distC, \@distN))))};
	
	$maxND=size(finalDist,2)-1;
	for (my $i = 3; $i < $maxND; $i++)
	{
		if ($finalDist[$i] < $obsCThreshold && $finalDist[$i-1] < $obsCThreshold && $finalDist[$i-2] >= $obsCThreshold)
		{
			$maxND = $i - 3;
			last;
		}
	}
	
	@distND = @finalDist(1:(maxND+1));
}

my $intype = ['Comma Seperated Values', '.csv'];
my $outtype = ['Text Files', ['.txt', '.text']];
my $infilename = "", $outfilename = "";

my $mw = new MainWindow;

$infilename = $mw->getOpenFile(-filetypes=>$intype);
if (-e $infilename && -r $infilename)
{
	print DEBUG "Input file is $infilename\n";
	open (my ifh, '<', $infilename) or die ('Could not open input file ' . $infilename . '.\n');
	$outfilename = $mw->getSaveFile(-filetypes=>$outtype, -initialfile=>"spif");
	if (-w $outfilename)
	{
		if ($outfilename eq "") {$outfilename = "spif.txt";}
		print DEBUG "Output file will be $outfilename\n";
		open (my ofh, '>', $outfilename) or die ('Could not open output file ' . $outfilename . '.\n');
		
		#do all the things to ifh
		while (<ifh>)
		{
			#do all the things to $_
			if ($_ =~ m/^\d.+$/)	#as long as $_ has data
			{
				my $start, $end, $z, $mz, $maxND, $maxD, $rt, $score, $delta, @distND, $seq, $title;
				#match variables to $_
				/(?:[^,]+,){6}(<mz>?[^,]+),(<z>?[^,]+),(?:[^,]+),(<delta>?[^,]+),(<start>?[^,]+),(<end>?[^,]+),(<score>?[^,]+),(?:[^,]+,)(<seq>?[^,]+),(?:[^,]+,)(<title>?[^,]+)/;
				$start = \g{start};
				$end = \g{end};
				$z = \g{z};
				$mz = \g{mz};
				$score = \g{score};
				$delta = \g{delta};
				$seq = \g{seq};
				$title = \g{title};
				$title =~ /^.+?RT:(<rt>?\d+(?:\.\d+))/;
				$rt = \g{rt};
				#call pepinfo
				&pepinfo($seq, \$maxND, \$maxD, \@distND);
				#print output to ofh
				print ofh, ($start . '\t' . $end . '\t' . $z . '\t' . $mz . '\t' . $maxND . '\t' . $maxD . '\t' . $rt . '\t' . $score . '\t' . $delta . '\t' . @distND . '\n');
			}
		}	
	}
	else
	{
		if (!-w $outfilename) {die ("Output file " . $outfilename . " is not writeable.\n";)}
	}
}
else
{
	if (-e $infilename) {die ("Input file " . $infilename . " does not exist.\n";)}
	elsif (-r $infilename) {die ("Input file " . $infilename . " is not readable.\n";)}
	else {die "This should be unreachable.\n"}
}
close DEBUG;