<?php
$option = $_GET['params'];
$defaultFile = $_GET['defaultFile'];
$newFile = $_GET['newFile'];


$optionArray = explode("!",$option);
$handle = fopen($newFile,"w+");
if(!$handle)
{
    echo("Could not open text of file");
    return;
}
#write the new stuff
$size = count($optionArray);
echo("SIZE=" . $size);
for($i = 0; $i < $size; $i++)
{
    echo($optionArray[$i]);
    $val = fwrite($handle,$optionArray[$i]);
    
    if(!$val)
    {
       echo("ERROR: Could not write to file");
      return;
    }
    fwrite($handle,"\n");
}
#copy over the old stuff
$oldVals = file_get_contents($defaultFile);

if(!$oldVals)
{
    echo("ERROR: Unable to open/read default params file");
    return;
}
$val = fwrite($handle,$oldVals);
if(!$val)
{
    echo("ERROR: Could not write to file");
    return;
}
$val = fclose($handle);
if(!$val)
{
    echo("ERROR: Unable to close file");
    return;
}
echo("Created new file " . $newFile);
?>