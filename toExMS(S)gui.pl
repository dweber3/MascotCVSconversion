#		toExMS(S)gui.pl
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