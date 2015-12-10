<?php
$option = $_GET['option'];
$file = $_GET['paramFile'];

#$handle = fopen($file,"r");
$handle = file_get_contents($file);
if(!$handle)
{
    echo("Could not read text of file");
    return;
}
$pos = strpos($handle,$option,0);
if(!$pos)
{
    echo("COULD NOT FIND GUY!!!!");
    return;
}
$pos = $pos + strlen($option) + 1;
$endPos = strpos($handle,"\n",$pos);

$val = substr($handle, $pos,$endPos-$pos);
echo($val);
?>