#!/usr/bin/perl

use CGI  qw(:standard);
use URI::Escape;

$| = 1;

$SPS_DIR = "";
$SPS_BIN = $SPS_DIR . "bin/";


$ENV{'LD_LIBRARY_PATH'} = "";

$cin = <> . "&" . $ENV{'QUERY_STRING'};

print("Content-type: text/html\r\n\r\n");


if(length($cin) > 0) {

  $call = $SPS_BIN . "spsReportServer.exe ";

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


  # Uncomment the following line to log the calls
  # `echo $call >> log_spsReports.txt`;


  # Below are a set of 3 lines, with different logging options
  # . Only one should be uncommented.

  # Uncomment the following code line to log errors
  #$cout = `$call 2> log_error_spsReports.txt`;

  # Uncomment the following code line to log specplot's output
  #cout = `$call  2>&1 1> log_all_spsReports.txt`;

  # Uncomment the following line to suppress all errors (default)
  $cout = `$call`  2>&1 1>/dev/null`;


  print($cout);
}
