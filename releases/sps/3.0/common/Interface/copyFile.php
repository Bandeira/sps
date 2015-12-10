<?php
$src = $_GET['src'];
$dest = $_GET['dest'];

$oldumask = umask(0);

if(!copy($src,$dest))
    echo("ERROR: Unable to copy file!");
umask($oldumask);
chmod($dest, 0777);
echo(realpath($dest));

?>