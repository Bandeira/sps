#!/usr/bin/perl


sub generate_random_string
{
	my $length_of_randomstring=shift;# the length of
			 # the random string to generate

	my @chars=('a'..'z','A'..'Z','0'..'9','_','-','+');
	my $random_string;
	foreach (1..$length_of_randomstring)
	{
		# rand @chars will generate a random
		# number between 0 and scalar @chars
		$random_string.=$chars[rand @chars];
	}
	return $random_string;
}



use CGI  qw(:standard);
use URI::Escape;

$| = 1;

$TMP = "/tmp";
$SPS_DIR = "";
$SPS_BIN = $SPS_DIR . "bin/";


$ENV{'LD_LIBRARY_PATH'} = "<>";

$cin = <> . "&" . $ENV{'QUERY_STRING'};

print("Content-type: image/png\r\n\r\n");

if(length($cin) > 0) {

  $call = $SPS_BIN . "specplot ";

  @params = split(/&/, $cin);
  foreach $param (@params) {
    ($key, $value) = split(/=/, $param);
    $key = uri_unescape($key);
    $value = uri_unescape($value);
    $call .= " " . $key;

    @values = split(/\+/, $value);
    foreach $value (@values) {
      #$value = quotemeta $value;
      $call .= " " . $value;
    }

  }

  $call .= " --outdir " . $TMP;

  $call .= " --outfile " . &generate_random_string(80);


  # Uncomment the following line to log the calls
  # `echo $call >> log_specplot.txt`;


  # Below are a set of 3 lines, with different logging options
  # . Only one should be uncommented.

  # Uncomment the following code line to log errors
  #$cout = `$call 2> log_error_specplot.txt`;

  # Uncomment the following code line to log specplot's output
  #cout = `$call  2>&1 1> log_all_specplot.txt`;

  # Uncomment the following line to suppress all errors (default)
  $cout = `$call`  2>&1 1>/dev/null`;


  print($cout);
}
