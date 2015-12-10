<?php
$fileName = $_GET['fileName'];
$seqName = $_GET['seqName'];
$sequence = $_GET['seq'];

$handle = fopen($fileName,"a");

if($handle == false)
  echo "2";

$val = fwrite($handle,">" . $seqName . "\n");
if($val == false)
  echo "3";

$val = fwrite($handle,$sequence . "\n");
if($val == false)
  echo "4";

$val = fclose($handle);
if($val == false)
  echo "5";

echo "1";
?>