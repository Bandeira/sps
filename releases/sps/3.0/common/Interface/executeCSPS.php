<?php
$paramFile = $_GET['paramFile'];
$exeDir = $_GET['exeDir'];
$projDir = $_GET['projectDir'];

# ini_set('max_execution_time', 900); //15 minutes

$file = $projDir . "/status.txt";
if (file_exists($file)){
  return;
 }
if(!chdir($projDir))
  echo($projDir);
$command = $exeDir . "/main_specnets " . $paramFile . " -ll 0 -lf log.txt -m -g > stdout.txt &";
$oldumask = umask(0);
$val = exec($command, $output, $retval);
# for($i = 0; $i < count($output); $i++)
#	echo($output[$i] . " ");

umask($oldumask);

echo($retval);

?>