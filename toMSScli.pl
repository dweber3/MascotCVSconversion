#		toMSS.pl
#	
#	Arranges information in Mascot-exported CVS files to that necessary for
#	Mass Spec Studio. No need to export to .xlsx afterwards! Made to work with Mascot as run by the IBBR as of 2014-07-18.
#	Mascot claims last update 2010-03-30.
#	
###############################################################################
#
#		Necessary Export Options:
#	-Protein Hit Information
#	-Peptide Match Information
#		=Experimental Charge
#		=Start
#		=End
#		=Sequence
#		=Query Title
#
#		Sample line follows:
#	prot_hit_num,	prot_acc,	pep_query,	pep_rank,	pep_isbold,	pep_isunique,	pep_exp_mz,	pep_exp_z,	pep_calc_mr,	pep_delta,	pep_start,	pep_end,	pep_score,	pep_res_before,	pep_seq,	pep_res_after,	pep_scan_title
#	1,				NLGT,		581,		1,			1,			1,				520.2807,	1,			519.2727,		0.0007,		130,		134,		31.09,		A,				SVVCL,		L,				File:06262014-Fab-MSMS-1.mzXML Scans:801 RT:4.6696min Charge:1+ Fragmentation:cid
#
###############################################################################
#
#	Output Format:
#	ID,Protein,Start,Stop,Sequence,z,m/z,RT
#
###############################################################################
use strict;
use warnings;

#open (DEBUG, ">log.txt");

my (@start_end, @AoR, @row, $RT);
my $ID = -1;
my $Prot = -1;
my $Start = -1;
my $Stop = -1;
my $Seq = -1;
my $z = -1;
my $mz = -1;
my $RTl = -1;
$RT = 0;

my $infilename = "";
my $outfilename = "";
#get $infilename
print "Please enter the input file:\t";
$infilename = <STDIN>;
chomp $infilename;
#print $DEBUG "Input file is $infilename.\n";
if (-e $infilename && -r $infilename)
{
	open (my $ifh, '<', $infilename) or die ('Could not open input file ' . $infilename . '.\n');
	#get $outputfilename
	print "Enter desired output filename. If none is entered, file will be 'out.csv'.\t";
	$outfilename = <STDIN>;
	chomp $outfilename;
	if ($outfilename eq "") {$outfilename = "out.csv";}
	#print $DEBUG "Output file will be $outfilename.\n";
	open (my $ofh, '>', $outfilename) or die ('Could not open output file ' . $outfilename . '.\n');
	if (-w $outfilename)
	{
		while (<$ifh>)
		{
			#parsing input file here
			
			#possible modes:
			#	0:	initial state
			#	1:	body of csv, following lines contain almost all of the necessary information
			
			if ($_ =~ /^prot_hit.*/ )
			{
				#header row
				#print DEBUG $_;
				@row = split(/,/, $_);
				for (my $i = 0; $i < scalar @row; $i++)
				{
					#ID fields
					if ($row[$i] eq "pep_query") {$ID = $i;}
					elsif ($row[$i] eq "prot_acc") {$Prot = $i;}
					elsif ($row[$i] eq "pep_start") {$Start = $i;}
					elsif ($row[$i] eq "pep_end") {$Stop = $i;}
					elsif ($row[$i] eq "pep_seq") {$Seq = $i;}
					elsif ($row[$i] eq "pep_exp_z") {$z = $i;}
					elsif ($row[$i] eq "pep_exp_mz") {$mz = $i;}
					elsif ($row[$i] eq "pep_scan_title\n") {$RTl = $i;}
					else {}
			
					#print DEBUG "$i: $row[$i]\n\tID: $ID; Prot: $Prot; Start: $Start; Stop: $Stop; Seq: $Seq; Z: $z; M/Z: $mz; RT: $RTl;\n"
				}
				push @AoR, "ID,Protein,Start,Stop,Sequence,m/z,z,RT, RT Variance, XIC m/z, XIC Adjustment, XIC Width, Number of Isotopic Peaks,Notes\n";
			}
			elsif ($_ =~ /^\d.*/)
			{
				#body section
				@row = split(/,/, $_);
				#print "$row[18]\n";
				if ($row[$RTl] =~ /.*RT:(\d+(?:\.\d+)).*/) {$RT = $1;}
				chomp($row[$RTl]);
				push @AoR, join (",", $row[$ID], $row[$Prot], $row[$Start], $row[$Stop], $row[$Seq], $row[$z], $row[$mz], $RT, "\n");		
			}
			else
			{
				#nothing to see here, move along.
			}
		}
		close $ifh;

		#Debug prints
		#for (my $i = 1; $i < $#AoR; $i++)
		#{
			#print "$AoR[$i]";
		#}

		#Print to output file.
		for (my $i = 0; $i < $#AoR; $i++)
		{
			print $ofh "$AoR[$i]";
		}
		
		close $ofh;
	}
}