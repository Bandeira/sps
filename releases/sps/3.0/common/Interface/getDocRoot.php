<?php

$docRoot = getenv("DOCUMENT_ROOT");
if(strlen($docRoot) == 0)
  {
    $localpath=getenv("SCRIPT_NAME");
    $absolutepath=realpath($localPath);
    // a fix for Windows slashes
    $absolutepath=str_replace("\\","/",$absolutepath);
    $docroot=substr($absolutepath,0,strpos($absolutepath,$localpath));

  }
echo($docRoot);

?>