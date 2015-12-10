<?php
$appendToMe = $_GET['appendToMe'];
$otherFile = $_GET['otherFile'];

$handle = fopen($appendToMe,"a");
$inHandle = fopen($otherFile,"r");
if($handle == false || $inHandle == false)
  echo "2";

$fileStuff = fread($inHandle,filesize($otherFile));
if($fileStuff == false)
  echo "3";
$val = fclose($inHandle);
if($val == false)
  echo "4";
$val = fwrite($handle,$fileStuff);
if($val == false)
  echo "5";

$val = fclose($handle);
if($val == false)
  echo "6";
echo "1";
?>