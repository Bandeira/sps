<?php
$fileName = $_GET['fileName'];

if(file_exists($fileName))
  echo("1");
else
  echo("0");
?>