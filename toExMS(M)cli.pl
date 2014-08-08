#		toExMS(M)cli.pl
#	
#	Arranges information in Mascot-exported CVS files to that necessary for
#	ExMS. Remember to export to .xlsx afterwards! Made to work with Mascot as run by the IBBR as of 2014-07-18.
#	Mascot claims last update 2010-03-30.
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
#open (my $DEBUG, '>', "log.txt");

my $infilename = "";
my $outfilename = "";
#get $infilename
print "Please enter the input file:\t";
$infilename = <STDIN>;
chomp $infilename;
#print $DEBUG "Input file is $infilename.\n";
if (-e $infilename && -r $infilename)
{
	open (my $ifh, '<', $infilename) or print $DEBUG "Could not open input file $infilename.\n";
	
	#get $outputfilename
	print "Enter desired output filename. If none is entered, file will be 'out.csv'.\t";
	$outfilename = <STDIN>;
	chomp $outfilename;
	if ($outfilename eq "") {$outfilename = "out.csv";}
	#print $DEBUG "Output file will be $outfilename.\n";
	if (open (my $ofh, '>', $outfilename))
	{
		if (!-w $outfilename) {die "Output file $outfilename is not writable.\n";}
		my (@start_end, @AoR, $RT, @row);
		my $mode = 0;
		my $hit = -1;
		my $prot = -1;
		my $mz = -1;
		my $z = -1;
		my $mr = -1;
		my $delta = -1;
		my $score = -1;
		my $start = -1;
		my $end = -1;
		my $bef = -1;
		my $seq = -1;
		my $aft = -1;
		my $scan = -1;
		$RT = 0;

		#possible modes:
		#	0:	initial state
		#	1:	body of csv, following lines contain almost all of the necessary information

		while (<$ifh>)
		{
			#parsing input file here
			if ($_ =~ /^prot_hit.*/ )
			{
				#header row
				$mode = 1;
				#print DEBUG $_;
				@row = split(/,/, $_);
				for (my $i = 0; $i < scalar @row; $i++)
				{
					#ID fields
					if ($row[$i] eq "prot_hit_num") {$hit = $i;}
					elsif ($row[$i] eq "prot_acc") {$prot = $i;}
					elsif ($row[$i] eq "pep_exp_mz") {$mz = $i;}
					elsif ($row[$i] eq "pep_exp_z") {$z = $i;}
					elsif ($row[$i] eq "pep_calc_mr") {$mr = $i;}
					elsif ($row[$i] eq "pep_delta") {$delta = $i;}
					elsif ($row[$i] eq "pep_score") {$score = $i;}
					elsif ($row[$i] eq "pep_start") {$start = $i;}
					elsif ($row[$i] eq "pep_end") {$end = $i;}
					elsif ($row[$i] eq "pep_res_before") {$bef = $i;}
					elsif ($row[$i] eq "pep_seq") {$seq = $i;}
					elsif ($row[$i] eq "pep_res_after") {$aft = $i;}
					elsif ($row[$i] eq "pep_scan_title\n") {$scan = $i;}
					else {}
				}
				push @AoR, "prot_hit_num,prot_acc,pep_exp_mz,pep_exp_z,pep_calc_mr,pep_delta,pep_score,retention time,pep_res_before,pep_seq,pep_res_after,pep_scan_title\n";
			}
			elsif ($_ =~ /^\d.*/ && $mode == 1)
			{
				#body section
				@row = split(/,/, $_);
				#print "$row[16]\n";
				if ($row[$scan] =~ /.*RT:(\d+(?:\.\d+)).*/) {$RT = $1;}
				chomp($row[$scan]);
				push @AoR, join ("", (join (',', ($row[$hit], $row[$prot], $row[$mz], "", $row[$z], $row[$mr], $row[$delta], $row[$score], $RT, $row[$bef], $row[$seq], $row[$aft], $row[$scan], " ", " ", " ", " " . $row[$start] . "-" . $row[$end])), "\n"));
			}
			else
			{
				#nothing to see here, move along.
				$mode = 0;
			}
		}
		close $ifh;

		#Debug prints
		#for (my $i = 1; $i < $#AoR; $i++)
		#{
			#print DEBUG "$AoR[$i]";
		#}

		#Print to output file.
		for (my $i = 0; $i < $#AoR; $i++)
		{
			print $ofh "$AoR[$i]";
		}

		close $ofh;
	}
	else {print $DEBUG "Could not open output file $outfilename.\n";}
}
else
{
	if (-e $infilename) {die "Input file $infilename does not exist.\n";}
	elsif (-r $infilename) {die "Input file $infilename is not readable.\n";}
	else {die "This should be unreachable.\n"}
}