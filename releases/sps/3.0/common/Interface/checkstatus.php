<?php
$paramFile = $_GET['paramFile'];
$exeDir = $_GET['exeDir'];
$projDir = $_GET['projectDir'];
$file = $projDir . "/status.txt";
if (file_exists($file)){
	readfile($file);
}else{
	echo("Error (PHP): Unable to determine status. File not found.");
	
}
?>