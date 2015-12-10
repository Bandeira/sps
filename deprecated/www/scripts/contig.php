<?php

$PATH=/home/pbouchard/bin

$N=28
$prefix="test"
$outdir="/data/sites/projects/sps/test"

print <<<EOF
<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN"
  "http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">

<HTML xmlns="http://www.w3.org/1999/xhtml">

<HEAD>

  <meta http-equiv="Content-Type" content="text/html; charset=ISO-8859-1" />
  <title>UCSD Computational Mass Spectrometry Website</title>
  <link href="../styles/main.css" rel="stylesheet" type="text/css" />
  <link rel="shortcut icon" href="../images/favicon.ico" type="image/icon" />
  <script src="../scripts/util.js" language="javascript" type="text/javascript"></script>
  <script src="../scripts/download.js" language="javascript" type="text/javascript"></script>
  <script src="../scripts/render.js" language="javascript" type="text/javascript"></script>
  <script src="../scripts/inspect.js" language="javascript" type="text/javascript"></script>

  <SCRIPT src="../sorttable.js"></script>
  <SCRIPT type="text/javascript">
  <!--hide

  function modalWin(filename)
  {
    if (window.showModalDialog)
    {
      window.showModalDialog(filename, "name", "dialogWidth:640px; dialogHeight:480px");
    }
    else
    {
      window.open(filename, 'name', 'width=640, height=480, toolbar=no, directories=no, status=no, menubar=no, scrollbars=no, resizable=no, modal=yes');
    }
  }

//-->
</SCRIPT>
</HEAD>
  <div id="rb_logos">
    <h3 align=center><img src="../pagelogo_ccms_left.jpg" border="0"/></h3>
  </div>
<BODY>
<p>This is a test.</p>
EOF;

system("$PATH/abruijn_test -s example2/stars.pklbin -c example2/component_info.bin -o example2/sps_seqs.pklbin -n $N --outdir $outdir --prefix $prefix --png --html")
system("$PATH/gnuplot $outdir/$prefix.$N.gnu")

system("$PATH/spectrum_test -s example2/specs_ms.pklbin -p $outdir/$prefix.$N.txt -z .3 --outdir $outdir --prefix $prefix --png --html")
system("$PATH/gnuplot $outdir/$prefix.$N.gnu")

system("$PATH/spectrum_test -s example2/specs_ms.pklbin -p $outdir/$prefix.$N.txt -z 1. --outdir $outdir --prefix l_${prefix} --png")
system("$PATH/gnuplot $outdir/l_${prefix}.$N.gnu")

print <<<EOF
</BODY>
EOF;

?>
