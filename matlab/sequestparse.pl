#! /usr/bin/perl
$count=0;
open (TMP,">TMP") or die;
#trim down the file
while (<>) {
#count total before cleaning up junk scans
if ($_ =~ /\d*  <A TARGET/) {
		$tt++;
}

$_ =~ s/<FONT COLOR=\"#DD00DD\">//g;
$_ =~ s/<A href=\"http:\/\/none\" target = Win1>//g;
$_ =~ s/<A TARGET="Win1" HREF=\"\/none\">//g;
$_ =~ s/<\/FONT>//g;
$_ =~ s/\*//g;
$_ =~ s/<\/FONT>//g;
$_ =~ s/<A href = \"http: \/\/none\" target = Win1>//g;
$_ =~ s/<\/A>//g;
$_ =~ s/<A TARGET=\"Win1\" HREF=\"\/cgi-shl\/web_showoutput.exe\?OutFile=C:\\Xcalibur\\results\/LCQ.*a\d*\.\d*\.\d*\.\d\.out\">//g;
$_ =~ s/<A TARGET=\"Win1\" HREF=\"\/cgi-shl\/web_display\.exe\?Dta=C\:\\Xcalibur\\results\/LCQ.*a\d*\.\d*\.\d*\.\d\.dta&amp;MassType=1&amp;NumAxis=1&amp;Pep=[A-Z]*\">//g;
$_ =~ s/<A TARGET=\"Win1\" href=\"http:\/\/www.ncbi.nlm.nih.gov\/blast\/blast.cgi.*isset\">//g;
$_ =~ s/<A TARGET.*web_retrieve.exe.*MassType=.\">//g;
$_ =~ s/<FONT COLOR=\"DD00DD\">//g;
$_ =~ s/<\/a>//g;
$_ =~ s/<HTML>//g;
$_ =~ s/<A TARGET.*BLAST.*>//g;


if ($_ =~ /spectrum/) {
	$_ ="";
}
if (defined $_) {
	print TMP "$_";
	$CC++;
}
}
close TMP;
print STDERR "Finished tmp write to\n";
open (TMPIN,"<TMP")or die;
while (<TMPIN>) {
	if ($_ =~ /\d*/) {
		$tl++;
	}
}
close TMPIN;

open (TMPIN,"<TMP") or die;
print "Line,file,scanstart,scanend,chargestate,mass,xcorr,deltacn,sp,Rsp,ions,reference,MW,sequence\n";
while (<TMPIN>) {																									     
	  																																			#ions--------------		gi----------------      MW--	sequence-------------------------     gi and protein name						
	   	if ( $_ =~ /(\d*)\s*([A-Za-z0-9\_]*LCQ\d*a\d{2,3})\.(\d*)\.(\d*)\.(\d)\s*(\d*\.\d*)\s*\(.\d*\.\d\)\s*(\d*\.\d*)\s*(\d*\.\d*)\s*(\d*\.\d*)\s*(\d*)\s*(\d*)(\/)\s{0,3}(\d*)\s*([A-Za-z0-9\|\_\.\+]*)\s*(\d*)\s*([A-Za-z0-9\-\.\|\[\]\s]*\.[A-Z\-])\s*([A-Za-z0-9\-\.\s\[\]\<\>\|]*)/){#([A-Za-z0-9\-\.\s\]\[\<\>\|]*)\s*/) {
			$bh++;
			$ln=$1;
			$file=$2;
			$scans=$3;
			$scanf=$4;
			$cs=$5;
			$mass=$6;
			$xcorr=$7;
			$deltacn=$8;
			$sp=$9;
			$Rsp=$10;
			$ions=$11;
			$ionslash=$12;
			$ione = $13;
			$ref=$14;
			$mw=$15;
			$seqs=$16;
			$seqm=$17;
			$seqe=$18;
			chomp($seqs,$seqm,$seqe,$ref);
			chop ($seqm);
			
			
			print "$ln,$file,$scans,$scanf,$cs,$mass,$xcorr,$deltacn,$sp,$Rsp,$ions\\$ione,$ref,$mw,$seqs,$seqm\n";
	   }
	   else {
	   	if ($_ !~ /^[<MLC\-\.\d]/ and $_ !~ /A TARGET/ and $_ !~ /\-*/) {
			push (@mis,$_);
			$mis++;
		}
	   }
	   
		
	
}
close TMPIN;

print STDERR "Total lines before cleaning $tt\n";
print STDERR "Total usable lines:$bh\n";
#print STDERR "# of lines copied:$CC\n";
#$mis=$tl-$bh;
if ($mis > 0) {
	print STDERR "Missed lines: $mis\n";
	print STDERR "Check MISSED file\n";
}
open (MIS,">MISSED") or die;
foreach $ms (@mis) {
	print MIS "$ms";
}
close MIS;
