#MascotCVSconversion
===================

##perltest.pl
###Warning: Permanently unstable!
Meant as an internal place to develop subroutines before copying them to scripts that actually put them to use. Always assume that it is less stable than property on top of a currently active volcano.

##Useful scripts:
Dependencies: | Required by: | Where to find them:
--- | --- | ---
Perl | All Perl scripts | http://www.perl.org/
Tk  | Only GUI versions | Acquired through the command _cpanm Tk_

##ExMS(M)     = ExMS (Matlab module)

Arranges information in Mascot-exported CVS files to that necessary for	ExMS. Remember, one needs to export it to .xlsx afterwards! Made to work with Mascot as run by the IBBR as of 2014-07-18. Mascot claims last update 2010-03-30.

####toExMS(M)cli.pl	= Command Line Interface version. Requires typing in filenames.
####toExMS(M)gui.pl	= Graphical User Interface version. Opens windows for the user to select files. Has additional dependencies.

####Required Export options:
	(If it isn't listed, this doesn't care about it.)

	Protein Hit Information

	Peptide Match Information
		Experimental Charge
		Calculated Mr (Da)
		Mass error (Da)
		Start
		End
		Score
		Sequence
		Query title

	Sample line follows:
	'''
	prot_hit_num,prot_acc,pep_query,pep_rank,pep_isbold,pep_isunique,pep_exp_mz,pep_exp_z,pep_calc_mr,pep_delta,pep_start,pep_end,pep_score,pep_res_before,pep_seq,pep_res_after,pep_scan_title
	1,NLGT,581,1,1,1,520.2807,1,519.2727,0.0007,130,134,31.09,A,SVVCL,L,File:06262014-Fab-MSMS-1.mzXML Scans:801 RT:4.6696min Charge:1+ Fragmentation:cid
	'''

##ExMS(S)  = ExMS (Standalone)

Arranges information in Mascot-exported CVS files to that necessary for ExMS. No need to export to .xlsx afterwards! Made to work with Mascot as run by the IBBR as of 2014-07-18. Mascot claims last update 2010-03-30.

####toExMS(S)cli.pl	= Command Line Interface version. Requires typing in filenames.
####toExMS(S)gui.pl	= Graphical User Interface version. Opens windows for the user to select files. Has additional dependencies.

####Required Export options:
	(If it isn't listed, this doesn't care about it.)

	Protein Hit Information

	Peptide Match Information
		Experimental Charge
		Calculated Mr (Da)
		Mass error (Da)
		Start
		End
		Score
		Sequence
		Query title

	Sample line follows:
	'''
	prot_hit_num,prot_acc,pep_query,pep_rank,pep_isbold,pep_isunique,pep_exp_mz,pep_exp_z,pep_calc_mr,pep_delta,pep_start,pep_end,pep_score,pep_res_before,pep_seq,pep_res_after,pep_scan_title
	1,NLGT,581,1,1,1,520.2807,1,519.2727,0.0007,130,134,31.09,A,SVVCL,L,File:06262014-Fab-MSMS-1.mzXML Scans:801 RT:4.6696min Charge:1+ Fragmentation:cid
	'''

##toMSS.pl       = Mass Spec Studio

Work in Progress. 

Dependencies: | Required by: | Where to find them:
--- | --- | ---
Perl | All Perl scripts | http://www.perl.org/
Tk  | Only GUI versions | Acquired through the command _cpanm Tk_

##MSS      = Mass Spec Studio

Arranges information in Mascot-exported CVS files to that necessary for	ExMS. Remember to export to .xlsx afterwards! Made to work with Mascot as run by the IBBR as of 2014-07-18. Mascot claims last update 2010-03-30.

####toMSScli.pl	= Command Line Interface version. Requires typing in filenames.
####toMSSgui.pl	= Graphical User Interface version. Opens windows for the user to select files. Has additional dependencies.

####Required Export options:
	(If it isn't listed, this doesn't care about it.)

	Protein Hit Information

	Peptide Match Information
		Experimental Charge
		Calculated Mr (Da)
		Mass error (Da)
		Start
		End
		Score
		Sequence
		Query title

	Sample line follows:
	'''
	prot_hit_num,prot_acc,pep_query,pep_rank,pep_isbold,pep_isunique,pep_exp_mz,pep_exp_z,pep_calc_mr,pep_delta,pep_start,pep_end,pep_score,pep_res_before,pep_seq,pep_res_after,pep_scan_title
	1,NLGT,581,1,1,1,520.2807,1,519.2727,0.0007,130,134,31.09,A,SVVCL,L,File:06262014-Fab-MSMS-1.mzXML Scans:801 RT:4.6696min Charge:1+ Fragmentation:cid
	'''

ExMS is available here: http://hx2.med.upenn.edu/download.html
Copyright 2009-2013 University of Pennsylvania

MSS is located here: http://structurems.ucalgary.ca/software/

