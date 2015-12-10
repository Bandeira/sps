<?php
$dirToMake = $_GET['dirToMake'];

$oldumask = umask(0);

if(!mkdir($dirToMake,0777,true))
    echo("ERROR: Unable to make directory!");
umask($oldumask);
chmod($dirToMake, 0777);
echo "1";
?>