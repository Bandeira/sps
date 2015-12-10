<?php
$dir = $_POST['dir'];
$root = "../data/"

if(file_exists ($root . $dir))
    echo "true";
else
    echo "false";

?>