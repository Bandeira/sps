
var browseRoot = "/";
var dataRoot = "/";


var defaultParamsFile="default.sps.params";
var newParamFileName = "genoms.sps.params";

var ar_ext = ['mzXML', 'mgf'];        // array with allowed extensions
var ar_ext_seq = ['fasta'];           //array with allowed sequence file extensions

var floatRegExp = /^[\+\-]?(\d+(\.\d*)?|\.?\d+)$/;
function getDocRoot()
{
    var xmlhttp;
    if (window.XMLHttpRequest)
    	{// code for IE7+, Firefox, Chrome, Opera, Safari
	    xmlhttp=new XMLHttpRequest();
    	}
    else
    	{// code for IE6, IE5
	    xmlhttp=new ActiveXObject("Microsoft.XMLHTTP");
    	}

    xmlhttp.open("Get","getDocRoot.php",false);
    xmlhttp.send(null);
	
    return xmlhttp.responseText;

}

//Utility functions

function fileExists(fileName)
{
    var xmlhttp;
    if (window.XMLHttpRequest)
    	{// code for IE7+, Firefox, Chrome, Opera, Safari
	    xmlhttp=new XMLHttpRequest();
    	}
    else
    	{// code for IE6, IE5
	    xmlhttp=new ActiveXObject("Microsoft.XMLHTTP");
    	}
    var params = "fileName=" + fileName;
    xmlhttp.open("Get","fileExists.php?" + params,false);
    xmlhttp.send(null);
    var resp = xmlhttp.responseText;
    
    if(resp == "1")
      return true;
    return false;
}

function makeDir(dirName)
{
     var xmlhttp;
    if (window.XMLHttpRequest)
    {// code for IE7+, Firefox, Chrome, Opera, Safari
	xmlhttp=new XMLHttpRequest();
    }
    else
    {// code for IE6, IE5
	xmlhttp=new ActiveXObject("Microsoft.XMLHTTP");
    }
    var params = "dirToMake=" + dirName;
    xmlhttp.open("Get","makeDir.php?" + params,false);
    xmlhttp.send(null);

if(xmlhttp.responseText == "1")
	return true;
else
	return false;
}



function getSeqFileName(chain,region)
{
	var fileName;
	           if(chain == "HC")
	           {
 			if(region == "V")
				fileName = chain + ".1";
			else if(region == "D")
				fileName = chain + ".2";
			else if(region == "J")
				fileName = chain + ".3";
			else if(region == "C")
				fileName = chain + ".4";
			else if(region == "W" || region == "Whole")
				fileName = chain + "_WholeSequence.fasta";
			else
			{
				alert("Invalid region '" + region + "!");
				return "";
			}
	    	   }
		   else if(chain == "LC")
		   {
			if(region == "V")
				fileName = chain + ".1";
			else if(region == "J")
				fileName = chain + ".2";
			else if(region == "C")
				fileName = chain + ".3";
			else if(region == "W" || region == "Whole")
				fileName = chain + "_WholeSequence.fasta";
			else
			{
				alert("Invalid region '" + region + "!");
				return "";
			}
		   }
		  else
		  {
			alert("Invalid chain '" + chain + "'!");
			return "";
		  }

	return fileName;
}

function writeSeqToFile(fileName,seqName,sequence)
{
    var xmlhttp;
    if (window.XMLHttpRequest)
    {// code for IE7+, Firefox, Chrome, Opera, Safari
	xmlhttp=new XMLHttpRequest();
    }
    else
    {// code for IE6, IE5
	xmlhttp=new ActiveXObject("Microsoft.XMLHTTP");
    }
    var params = "fileName=" + fileName + "&seqName=" + seqName + "&seq=" + sequence;
    xmlhttp.open("Get","addSequenceToFile.php?" + params,false);
    xmlhttp.send(null);

if(xmlhttp.responseText == "1")
	return true;
else
{
	alert(xmlhttp.responseText);
	return false;
	}

}

function appendToFile(appendToMe, otherFile)
{
    var xmlhttp;
    if (window.XMLHttpRequest)
    {// code for IE7+, Firefox, Chrome, Opera, Safari
	xmlhttp=new XMLHttpRequest();
    }
    else
    {// code for IE6, IE5
	xmlhttp=new ActiveXObject("Microsoft.XMLHTTP");
    }
    var params = "appendToMe=" + appendToMe + "&otherFile=" + otherFile;
    xmlhttp.open("Get","appendFiletoFile.php?" + params,false);
    xmlhttp.send(null);

    if(xmlhttp.responseText == "1")
	return true;
    else
{
	alert(xmlhttp.responseText);
	return false;
}
}


function basename(P)
{
    return P.replace(/\\/g,'/').replace( /.*\//, '' );
}

function ValidateForm(doc)
{
   //Do we have at least 1 spectrum file?
    var theTable = doc.getElementById("currentSpectrumTable");
    if(theTable.rows.length <= 1)
    {
	alert("Must specify at least 1 spectrum file!");
	return "null";
    }
    //Do we have at least one database checked?
    var HCCheckbox = doc.getElementById("HC_DB");
    var LCCheckbox = doc.getElementById("LC_DB");
        
    //Do we have a project directory
    projectString = doc.getElementById("pdir").value;
    if( projectString =="")
    {
	 alert("Please enter a project name");
	 return "null";
    }
    if(projectString.indexOf("..",0) >= 0 ||
       projectString.indexOf(" ",0) >= 0 ||
       projectString.indexOf("\\",0) >= 0 ||
       projectString.indexOf("/",0) >= 0)
    {
	alert("Invalid project name! Cannot contain '..', ' ', '/','\'");
	return "null";
    }
    
    //ValensProjectName = projectString;
   // alert(ValensProjectName);
    projectString = dataRoot + projectString;

    var projExists = fileExists(projectString);
    if(projExists)
    {
	alert("Project already exists, please choose another name");
	return "null";
    }
    return projectString;
}

function LaunchCSPS(doc)
{
    //Verify that we have all of the right info
    var projectString = ValidateForm(doc);
    if(projectString == "null")
      return;
 
    var xmlhttp;
    if (window.XMLHttpRequest)
    	{// code for IE7+, Firefox, Chrome, Opera, Safari
	    xmlhttp=new XMLHttpRequest();
    	}
    else
    	{// code for IE6, IE5
	    xmlhttp=new ActiveXObject("Microsoft.XMLHTTP");
    	}
	
    var params;
    
	
	
	if(!makeDir(projectString))
	{	
		alert("ERROR: Unable to make directory " + projectString);
		return;
	}	
	if(!makeDir(projectString + "/spectra"))
	{
		alert("ERROR: Unable to make directory " + projectString + "/spectra");
		return;
	}	
    //copy the spectra

    //alert("Made project directories");
    var i = 0;
    var fileName;
    var parentTol = 0.0;
    var fragTol = 0.0;
    var paramsFile="INPUT_SPECS_MS=";
    
    var theTable = doc.getElementById("currentSpectrumTable");
    for(i = 1; i < theTable.rows.length; ++i)
    {
	fileName = theTable.rows[i].cells[0].innerHTML;
	currParentTol = parseFloat(theTable.rows[i].cells[1].innerHTML);
	currFragTol = parseFloat(theTable.rows[i].cells[2].innerHTML);
	if(currParentTol > parentTol)
	    parentTol = currParentTol;
	if(currFragTol > fragTol)
	    fragTol = currFragTol;
	params = "src=" + fileName +"&dest=" + projectString + "/spectra/" + basename(fileName);
	xmlhttp.open("Get","copyFile.php?" + params,false);
	xmlhttp.send(null);
	absFileName=xmlhttp.responseText;
	if(i > 1)
	    paramsFile += ";" + absFileName;
	else
	    paramsFile += absFileName;
    }


    //Since the new lines will be wiped outin the GET message, use a different delim
    var delim = "!";
    var baseProjectName = projectString;
    params = "file=" + projectString;
    xmlhttp.open("Get","getAbsPath.php?" + params,false);
    xmlhttp.send(null);
    projectString = xmlhttp.responseText;
    
    params = "option=EXE_DIR&paramFile=" + defaultParamsFile;
    xmlhttp.open("Get","getParamsFileValue.php?" + params,false);
    xmlhttp.send(null);
    exeDir = xmlhttp.responseText;
    dbDir = exeDir + "/DBs_GenoMS/";
    
    paramsFile += delim + "PROJECT_DIR=" + projectString;
    paramsFile += delim + "FIXEDMOD=" + doc.getElementById("cysteine").value;
    paramsFile += delim + "RELATIVE_DIR=0";
    paramsFile += delim + "REPORT_JOB=" + baseProjectName;
    paramsFile += delim + "REPORT_DIR_SERVER=" + projectString;
    paramsFile += delim + "REPORT_USER=anonymous";
    paramsFile += delim + "REPORT_DIR=" + projectString + "/report";
    paramsFile += delim + "TOLERANCE_PEAK=" + fragTol;
    paramsFile += delim + "TOLERANCE_PM=" + parentTol;
    
//var ContaminantsCheckBox = doc.getElementById("cCheck");
//if(ContaminantsCheckBox.checked)
//	paramsFile += delim + "CONTAMINANTS=" + dbDir + "contaminants.RS.trie";

    //Do we have at least one database checked?
    var HCCheckbox = doc.getElementById("HC_DB");
    var LCCheckbox = doc.getElementById("LC_DB");

    var seqStr = document.getElementById("seqList").innerHTML;
    var fileStr = document.getElementById("dbFileList").innerHTML;    

//If no check boxes are selected and no user seqs are added, then wtf?
if(!HCCheckbox.checked && !LCCheckbox.checked && (seqStr == null | seqStr == "") && (fileStr == null || fileStr == ""))
    {
	alert("Must specify one or more chains, or add user sequence!");
	return "null";
    }


//There are no additional user databases, so just use the defaults
    if((seqStr == null || seqStr == "") &&
	(fileStr == null || fileStr == ""))
	{
	var fastaFileName = "";
	var constraintFileName = "";
	var wholeSeqFileName = "";
    	if(HCCheckbox.checked && LCCheckbox.checked)
    	{
		fastaFileName = dbDir + "IMGT_20120213_HC_LC.fasta"; 
		constraintFileName = dbDir + "IMGT_20120213_HC_LC.constraints";
		wholeSeqFileName = dbDir + "Mouse_Rat_WholeSequence.fasta";
    	}
    	else if(HCCheckbox.checked)
    	{
		fastaFileName = dbDir + "IMGT_20120213_HC.fasta"; 
		constraintFileName = dbDir + "IMGT_20120213_HC.constraints";
		wholeSeqFileName = dbDir + "Mouse_Rat_HC_WholeSequence.fasta";

    	}
    	else if(LCCheckbox.checked)
    	{
		fastaFileName = dbDir + "IMGT_20120213_LC.fasta"; 
		constraintFileName = dbDir + "IMGT_20120213_LC.constraints";
		wholeSeqFileName = dbDir + "Mouse_Rat_LC_WholeSequence.fasta";

    	}

	params = "src=" + fastaFileName +"&dest=" + projectString + "/" + basename(fastaFileName);
	xmlhttp.open("Get","copyFile.php?" + params,false);
	xmlhttp.send(null);
	absFileName=xmlhttp.responseText;
	paramsFile += delim + "DBCOMBINED=" + absFileName;

	params = "src=" + constraintFileName +"&dest=" + projectString + "/" + basename(constraintFileName);
	xmlhttp.open("Get","copyFile.php?" + params,false);
	xmlhttp.send(null);
	absFileName=xmlhttp.responseText;
	paramsFile += delim + "TEMPLATECONSTRAINTFILE=" + absFileName;


	params = "src=" + wholeSeqFileName +"&dest=" + projectString + "/" + basename(wholeSeqFileName);
	xmlhttp.open("Get","copyFile.php?" + params,false);
	xmlhttp.send(null);
	absFileName=xmlhttp.responseText;
	paramsFile += delim + "FASTA_DATABASE=" + absFileName;

	

    	paramsFile += "\n";
    }
    else
	{
		
	   var hcArray = new Array(5);
	   hcArray[0] = false;
	   hcArray[1] = false;
    	   hcArray[2] = false;
	   hcArray[3] = false;
	   hcArray[4] = false;

	   var lcArray = new Array(4);
	   lcArray[0] = false;
	   lcArray[1] = false;
    	   lcArray[2] = false;
	   lcArray[3] = false;
	
           //First deal with the sequences
	   if(seqStr != null && seqStr != "")
	   {
		var seqArray = seqStr.split(";");
		var i;
	  	for(i = 0; i < seqArray.length; ++i)
		{
		   var currSeqArray = seqArray[i].split(":");
		   var chain = currSeqArray[1];
		   var region = currSeqArray[2];
		   if(chain == "HC")
		   {
		 	if(region == "V")
				hcArray[0] = true;
			else if(region == "D")
				hcArray[1] = true;
			else if(region == "J")
				hcArray[2] = true;
			else if(region == "C")
				hcArray[3] = true;
			else if(region == "Whole")
				hcArray[4] = true;
		   }		
		   else if(chain == "LC")
		   {
		 	if(region == "V")
				lcArray[0] = true;
			else if(region == "J")
				lcArray[1] = true;
			else if(region == "C")
				lcArray[2] = true;
			else if(region == "Whole")
				lcArray[3] = true;
		   }		
	   
                   var fileName = getSeqFileName(chain,region);
		   if(fileName == "")
			return;
		   fileName = projectString + "/" + fileName;
		   if(writeSeqToFile(fileName,currSeqArray[0],currSeqArray[3]) == false)
		   {
			alert("Unable to create sequence file '" + fileName + "'!");
			return;	
		   }
		}
	   }
	   if(fileStr != null && fileStr != "")
	   {  
		var seqArray = fileStr.split(";");
		var i;
	  	for(i = 0; i < seqArray.length; ++i)
		{
		   var currSeqArray = seqArray[i].split(":");
		   var chain = currSeqArray[1];
		   var region = currSeqArray[2];
		   if(chain == "HC")
		   {
		 	if(region == "V")
				hcArray[0] = true;
			else if(region == "D")
				hcArray[1] = true;
			else if(region == "J")
				hcArray[2] = true;
			else if(region == "C")
				hcArray[3] = true;
			else if(region == "Whole")
				hcArray[4] = true;
		   }		
		   else if(chain == "LC")
		   {
		 	if(region == "V")
				lcArray[0] = true;
			else if(region == "J")
				lcArray[1] = true;
			else if(region == "C")
				lcArray[2] = true;
			else if(region == "Whole")
				lcArray[3] = true;
		   }		
                   var fileName = getSeqFileName(chain,region);
		   if(fileName == "")
			return;
		   fileName = projectString + "/" + fileName;
		   if(!appendToFile(fileName,currSeqArray[0]))	
		   {
			alert("Unable to append '" + currSeqArray[0] + "' to sequence file '" + fileName + "'!");
			return;	
		   }
		}
	   }
	var weSearchHC = false;
	var weSearchLC = false;
	var j;
	for(j = 0; j < hcArray.length; ++j)
	{
		if(hcArray[j])
			weSearchHC = true;
	}
	for(j = 0; j < lcArray.length; ++j)
	{
		if(lcArray[j])
			weSearchLC = true;
	}

	//Now we have to check that we have sequence from every region of any chain that is used!
	if((!hcArray[0] && weSearchHC) || HCCheckbox.checked) 
	{
		var fileName = getSeqFileName("HC","V");
		fileName = projectString + "/" + fileName;
		if(!appendToFile(fileName,dbDir + "IMGT_20120213_HC.1"))
		{
			alert("Unable to append '" + dbDir + "IMGT_20120213_HC.1' to '" + fileName + "'!");
			return;
		}
	}
	
	if((!hcArray[1] && weSearchHC) || HCCheckbox.checked)  
	{
		var fileName = getSeqFileName("HC","D");
		fileName = projectString + "/" + fileName;
		if(!appendToFile(fileName,dbDir + "IMGT_20120213_HC.2"))
		{
			alert("Unable to append '" + dbDir + "IMGT_20120213_HC.2' to '" + fileName + "'!");
			return;
		}
	}
	
	if((!hcArray[2] && weSearchHC) || HCCheckbox.checked)  
	{
		var fileName = getSeqFileName("HC","J");
		fileName = projectString + "/" + fileName;
		if(!appendToFile(fileName,dbDir + "IMGT_20120213_HC.3"))
		{
			alert("Unable to append '" + dbDir + "IMGT_20120213_HC.3' to '" + fileName + "'!");
			return;
		}
	}

	if((!hcArray[3] && weSearchHC)|| HCCheckbox.checked)   
	{
		var fileName = getSeqFileName("HC","C");
		fileName = projectString + "/" + fileName;
		if(!appendToFile(fileName,dbDir + "IMGT_20120213_HC.4"))
		{
			alert("Unable to append '" + dbDir + "IMGT_20120213_HC.4' to '" + fileName + "'!");
			return;
		}
	}

	if((!hcArray[4] && weSearchHC)|| HCCheckbox.checked)    
	{
		var fileName = getSeqFileName("HC","W");
		fileName = projectString + "/" + fileName;
		if(!appendToFile(fileName,dbDir + "Mouse_Rat_HC_WholeSequence.fasta"))
		{
			alert("Unable to append '" + dbDir + "Mouse_Rat_HC_WholeSequence.fasta' to '" + fileName + "'!");
			return;
		}
	}

	if((!lcArray[0] && weSearchLC) || LCCheckbox.checked)   
	{
		var fileName = getSeqFileName("LC","V");
		fileName = projectString + "/" + fileName;
		if(!appendToFile(fileName,dbDir + "IMGT_20120213_LC.1"))
		{
			alert("Unable to append '" + dbDir + "IMGT_20120213_LC.1' to '" + fileName + "'!");
			return;
		}
	}

	if((!lcArray[1] && weSearchLC)  || LCCheckbox.checked)   
	{
		var fileName = getSeqFileName("LC","J");
		fileName = projectString + "/" + fileName;
		if(!appendToFile(fileName,dbDir + "IMGT_20120213_LC.2"))
		{
			alert("Unable to append '" + dbDir + "IMGT_20120213_LC.2' to '" + fileName + "'!");
			return;
		}
	}

	if((!lcArray[2] && weSearchLC)  || LCCheckbox.checked)   
	{
		var fileName = getSeqFileName("LC","C");
		fileName = projectString + "/" + fileName;
		if(!appendToFile(fileName,dbDir + "IMGT_20120213_LC.3"))
		{
			alert("Unable to append '" + dbDir + "IMGT_20120213_LC.3' to '" + fileName + "'!");
			return;
		}
	}

	if((!lcArray[3] && weSearchLC)  || LCCheckbox.checked)   
	{
		var fileName = getSeqFileName("LC","W");
		fileName = projectString + "/" + fileName;
		if(!appendToFile(fileName,dbDir + "Mouse_Rat_LC_WholeSequence.fasta"))
		{
			alert("Unable to append '" + dbDir + "Mouse_Rat_LC_WholeSequence.fasta' to '" + fileName + "'!");
			return;
		}
	}
	
	if(HCCheckbox.checked)
		weSearchHC = true;
	if(LCCheckbox.checked)
		weSearchLC = true;
	//Now add all of this to the params file
	if(weSearchHC && !weSearchLC)
	{	paramsFile += delim + "DBROOTNAME=" + projectString + "/HC";
		paramsFile += delim + "FASTA_DATABASE=" + projectString + "/HC_WholeSequence.fasta";
	}
	else if(!weSearchHC && weSearchLC)
	{
		paramsFile += delim + "DBROOTNAME=" + projectString + "/LC";
		paramsFile += delim + "FASTA_DATABASE=" + projectString + "/LC_WholeSequence.fasta";
	}
	else
	{
		paramsFile += delim + "DBROOTNAME=" + projectString + "/HC;" + projectString + "/LC";

		if(!appendToFile(projectString + "/HC_LC_WholeSequence.fasta",projectString + "/HC_WholeSequence.fasta"))
		{
			alert("Unable to append 'HC_WholeSequence.fasta' to 'HC_LC_WholeSequence.fasta'");
			return;
		}	
		if(!appendToFile(projectString + "/HC_LC_WholeSequence.fasta",projectString + "/LC_WholeSequence.fasta"))
		{
			alert("Unable to append 'LC_WholeSequence.fasta' to 'HC_LC_WholeSequence.fasta'");
			return;
		}	
		paramsFile += delim + "FASTA_DATABASE=" + projectString + "/HC_LC_WholeSequence.fasta";
	}
    
    }

    var newParamsFile = projectString + "/" + newParamFileName;
    //Create the input file
    params = "params=" + paramsFile + "&defaultFile=" + defaultParamsFile + "&newFile=" + newParamsFile;
    xmlhttp.open("Get","extendParamsFileValue.php?" + params,false);
    
    xmlhttp.send(null);
    
    LaunchValens(projectString,exeDir);
    
}




function LaunchValens(projectString, exeDir)
{

    var params = "paramFile=" + newParamFileName + "&exeDir=" + exeDir  + "&projectDir=" + projectString;
    var launchWindow = window.location = "launcher.php?" + params;
    //var launchWindow = window.open("launcher.php?" + params, "Valens Jobs",
//			"height=520,width=400,toolbar=0,location=0,directories=0," +
//			"status=1,menubar=0,scrollbars=yes,resizeable=0");
    
    //alert(launchWindow.document.readyState);
		//$(launchWindow.document).ready(function() {
  	// Handler for .ready() called.
    //});
    /*var xmlhttp;
    	if (window.XMLHttpRequest)
    	{// code for IE7+, Firefox, Chrome, Opera, Safari
			xmlhttp=new XMLHttpRequest();
    	}
    	else
    	{// code for IE6, IE5
			xmlhttp=new ActiveXObject("Microsoft.XMLHTTP");
    	}

		launchWindow.document.getElementById("status").rows[1].cells[0].innerHTML = basename(projectString);
		launchWindow.document.getElementById("status").rows[1].cells[1].innerHTML = Date();
		launchWindow.document.getElementById("status").rows[1].cells[2].innerHTML = "Running";
		launchWindow.focus();
		var params = "paramFile=" + newParamFileName + "&exeDir=" + exeDir  + "&projectDir=" + projectString;
   	//launch
    	xmlhttp.open("Get","executeCSPS.php?" + params,false);
    	xmlhttp.send(null);
    	
    	launchWindow.document.getElementById("status").rows[1].cells[2].innerHTML="Finished with Status:" + xmlhttp.responseText;
	launchWindow.focus();*/
}
   
function checkFileExt(el) 
{
	
	
	// - www.marplo.net
	// get the file name and split it to separe the extension

	var ar_name = el.split('.');

	// for IE - separe dir path from name
	//var ar_nm = ar_name[0].split('\\');
	//alert(ar_name[0] + "->" + ar_nm);
	//for(var i=0; i<ar_nm.length; i++)
	//	var nm = ar_nm[i];

	// add the name in 'to'
	//document.getElementById(to).value = nm;

	// check the file extension
	var re = 0;
	for(var i=0; i<ar_ext.length; i++) {
		if(ar_ext[i] == ar_name[ar_name.length-1]) {
			re = 1;
			break;
		}
	}

	// if re is 1, the extension is in the allowed list
	if(re==1) {
	// enable submit
	    return true;
	}
	return false;
		  
}
function checkFileExtSeq(el) 
{
	// - www.marplo.net
	// get the file name and split it to separe the extension

	var ar_name = el.split('.');

	// for IE - separe dir path from name
	//var ar_nm = ar_name[0].split('\\');
	//alert(ar_name[0] + "->" + ar_nm);
	//for(var i=0; i<ar_nm.length; i++)
	//	var nm = ar_nm[i];

	// add the name in 'to'
	//document.getElementById(to).value = nm;

	// check the file extension
	var re = 0;
	for(var i=0; i<ar_ext_seq.length; i++) {
		if(ar_ext_seq[i] == ar_name[ar_name.length-1]) {
			re = 1;
			break;
		}
	}

	// if re is 1, the extension is in the allowed list
	if(re==1) {
	// enable submit
	    return true;
	}
	return false;
		  
}
	
function checkFileNamesinTable(table, val, col)
{
    var i;
if(table == null || table.rows == null)
	return false;
    for(i = 1; i < table.rows.length; ++i)
	{
	    if(val == table.rows[i].cells[col].childNodes[0].data)
		return true;

	}
    return false;


}
	
function addSpecFileRow(doc) {
	
    var fileName = doc.getElementById("chosenSpecFile").innerHTML;
    //var description = doc.forms[0].desc.value;
	var parentMass = doc.forms[0].pmass.value;
	var ionTol = doc.forms[0].itol.value;
	
	var protease = doc.forms[0].protease.value;
	var cGroup = doc.forms[0].cysteine.value;
	
	
	//Validate the values

	if(fileName == null || fileName == "")
	{
		alert("Please enter a spectrum file");
		//form.filess.focus();
		return;
	}

	if(!checkFileExt(fileName))
	    {
		alert("Invalid spectrum file format.  Must be '.mgf' or '.mzXML'");
		return;
	    }


	//if(description == null || description == "")
	//{
	//	alert("Please enter a file description");
	//	//form.desc.focus();
	//	return;
	//}
	

	if (!floatRegExp.test(parentMass)) {
		alert("Parent mass tolerance must be a real number.");
		//form.pmass.focus();
		return;
	}
	var fPMTol = parseFloat(parentMass);
	
	//alert("checked stuffA");
	if(fPMTol > 2.5 || fPMTol < 0)
	{
		alert("Invalid parent mass tolerance.  Must be between 0 and 2.5 Da");
		//form.pmass.focus();
		return;
	}
	if (!floatRegExp.test(ionTol)) {
		alert("Ion tolerance must be a real number.");
		//form.itol.focus();
		return;
	}
	//alert("checked stuffB");
	var fITol = parseFloat(ionTol);
	if(fITol > 1.0 || fITol < 0)
	{
		alert("Invalid ion mass tolerance.  Must be between 0 and 1 Da");
		//form.itol.focus();
		return;
	}
	
	
	

	//alert("building table!!");
	
	var theTable = doc.getElementById("currentSpectrumTable");
	
	//Check that the file doesn't already exist in the table!!

	if(checkFileNamesinTable(theTable,fileName,0))
	    {
		alert("Spectrum file already exists in the table!");
		return;
	    }
	//theTable.setAttribute("border","1");
	//theTable.setAttribute("class","ptm");
	
	//If the table is empty then it only has 1 lines
	//1. 'No spectrum files have been loaded'
	/*if ( theTable.rows.length == 1)
	{
		//delete No spectrum files have been added
	    theTable.deleteRow(0);
	    var hRow = theTable.insertRow(0);
		
		//insert a new header
		var newCell
					
		// an inserted row has no cells, so insert the cells
		newCell = hRow.insertCell(0)
		// give this cell its own id
		newCell.id = newCell.uniqueID
		// display the row's id as the cell text
		newCell.innerHTML = "File Name"
		newCell.width = "269"
		newCell.align ="left"
		newCell.bgcolor ="#0000FF"
					
		// reuse cell var for second cell insertion
		newCell = hRow.insertCell(1)
		newCell.id = newCell.uniqueID
		newCell.innerHTML = "PM Tolerance"
		newCell.width = "145"
					
		// reuse cell var for third cell insertion
		newCell = hRow.insertCell(2)
		newCell.id = newCell.uniqueID
		newCell.innerHTML = "Ion tolerance"
		newCell.width = "145"
					
		// reuse cell var for fourth cell insertion
		newCell = hRow.insertCell(3)
		newCell.id = newCell.uniqueID
		newCell.innerHTML = "Protease"
		newCell.width = "145"
					
		// reuse cell var for fifth cell insertion
		newCell = hRow.insertCell(4)
		newCell.width = "50"
					
		newCell = hRow.insertCell(5)
		newCell.id = newCell.uniqueID
		    newCell.innerHTML="Trash";
		newCell.style.display ="none";
		newCell.width = "145"
		
		}*/
	var row = theTable.insertRow(theTable.rows.length);
	var fileCell = row.insertCell(0);
	var pTolCell = row.insertCell(1);
	var iTolCell = row.insertCell(2);
	var protCell = row.insertCell(3);
	var tempCell = row.insertCell(4);
	var ctrlCell = row.insertCell(5);
	
	var ctrlButton = document.createElement("img");
	ctrlButton.src = "images/recycle.jpg";
	ctrlButton.className = "selectable"
	ctrlButton.onclick = function() {
		theTable.deleteRow(row.rowIndex);
		return false;
	}


	

	ctrlCell.appendChild(ctrlButton);
	fileCell.innerHTML = fileName;
	pTolCell.innerHTML = parentMass;
	iTolCell.innerHTML = ionTol;
	protCell.innerHTML = protease;
	
	//initialize values
	doc.getElementById("chosenSpecFile").innerHTML="No spectrum file selected"
	//doc.all.cysteine.value="C,+57";
	//doc.all.protease.value="Trypsin";
	//doc.all.pmass.value="2";
	//doc.all.itol.value="0.5";
		
	
}
function addDBFileRow(doc) {
	
    var fileName = doc.getElementById("chosenDBFile").innerHTML;
    		
    var chain = determineFileChain(doc);
    var region = determineFileRegion(doc);

	//Validate the values

	if(fileName == null || fileName == "")
	{
		alert("Please enter a sequence file");
		//form.filess.focus();
		return;
	}

	if(!checkFileExtSeq(fileName))
	    {
		alert("Invalid spectrum file format.  Must be '.fasta'");
		return;
	    
	}
	if(chain == null || chain == "")
	{
		alert("Please select a chain!");
		return;
	}
	if(region == null || region == "")
	{
		alert("Please select a region!");
		return;
	}
	if(chain == "LC" && region == "D")
	{
		alert("Invalid region for light chain");
		return;
	}


	//alert("building table!!");
	
	var theTable = doc.getElementById("currentSeqFiles");
	
	//Check that the file doesn't already exist in the table!!

	if(checkFileNamesinTable(theTable,fileName,0))
	    {
		alert("Sequence file already exists in the table!");
		return;
	    }
	//theTable.setAttribute("border","1");
	//theTable.setAttribute("class","ptm");
	
	//If the table is empty then it only has 1 lines
	//1. 'No spectrum files have been loaded'

	var row = theTable.insertRow(theTable.rows.length);
	var fileCell = row.insertCell(0);
	var chainCell = row.insertCell(1);
	var regionCell = row.insertCell(2);
	var tempCell = row.insertCell(3);
	var ctrlCell = row.insertCell(4);
	
	
	var ctrlButton = document.createElement("img");
	ctrlButton.src = "images/recycle.jpg";
	ctrlButton.className = "selectable"
	ctrlButton.onclick = function() {
		theTable.deleteRow(row.rowIndex);
		return false;
	}


	

	ctrlCell.appendChild(ctrlButton);
	fileCell.innerHTML = fileName;
	chainCell.innerHTML = chain;
	regionCell.innerHTML = region;

	
	//initialize values
	doc.getElementById("chosenDBFile").innerHTML="No sequence file selected"
	doc.getElementById("chainGroupFH").checked = true;
	doc.getElementById("regionGroupFV").checked = true;	
	
	//doc.all.cysteine.value="C,+57";
	//doc.all.protease.value="Trypsin";
	//doc.all.pmass.value="2";
	//doc.all.itol.value="0.5";
		
	
}

function determineSeqChain(doc)
{
	if(doc.getElementById("chainGroupSH").checked)
		return "HC";
	else if(doc.getElementById("chainGroupSL").checked)
		return "LC";
	return "";

}

function determineSeqRegion(doc)
{
	if(doc.getElementById("regionGroupSV").checked)
		return "V";
	else if(doc.getElementById("regionGroupSD").checked)
		return "D";
	else if(doc.getElementById("regionGroupSJ").checked)
		return "J";
	else if(doc.getElementById("regionGroupSC").checked)
		return "C";
	else if(doc.getElementById("regionGroupSW").checked)
		return "Whole";
	return "";

}


function determineFileChain(doc)
{
	if(doc.getElementById("chainGroupFH").checked)
		return "HC";
	else if(doc.getElementById("chainGroupFL").checked)
		return "LC";
	return "";

}

function determineFileRegion(doc)
{
	if(doc.getElementById("regionGroupFV").checked)
		return "V";
	else if(doc.getElementById("regionGroupFD").checked)
		return "D";
	else if(doc.getElementById("regionGroupFJ").checked)
		return "J";
	else if(doc.getElementById("regionGroupFC").checked)
		return "C";
	else if(doc.getElementById("regionGroupFW").checked)
		return "Whole";
	return "";

}



// return the value of the radio button that is checked
// return an empty string if none are checked, or
// there are no radio buttons

function addDBSeqRow(doc) {
	
    var seqName = doc.getElementById("sName").value;
    var seq = doc.getElementById("sVal").value;

	
    var chain = determineSeqChain(doc);
    var region = determineSeqRegion(doc);

	//Validate the values
	


	if(seqName == null || seqName == "")
	{
		alert("Please enter a sequence name!");
		//form.filess.focus();
		return;
	}

	if(seq == null || seq == "")
	{
		alert("Please enter a sequence!");
		//form.filess.focus();
		return;
	}

	if(chain == null || chain == "")
	{
		alert("Please select a chain!");
		return;
	}
	if(region == null || region == "")
	{
		alert("Please select a region!");
		return;
	}
	if(chain == "LC" && region == "D")
	{
		alert("Invalid region for light chain");
		return;
	}
	var theTable = doc.getElementById("currentSeqs");
	
	//Check that the file doesn't already exist in the table!!

	if(checkFileNamesinTable(theTable,seqName,0))
	    {
		alert("A sequence by that name already exists in the table!");
		return;
	    }
	
	//If the table is empty then it only has 1 lines
	//1. 'No spectrum files have been loaded'

	var row = theTable.insertRow(theTable.rows.length);
	var nameCell = row.insertCell(0);
	var chainCell = row.insertCell(1);
	var regionCell = row.insertCell(2);
	var tempCell = row.insertCell(3);
	var ctrlCell = row.insertCell(4);
	var hiddenSeqCell = row.insertCell(5);
	
	
	var ctrlButton = document.createElement("img");
	ctrlButton.src = "images/recycle.jpg";
	ctrlButton.className = "selectable"
	ctrlButton.onclick = function() {
		theTable.deleteRow(row.rowIndex);
		return false;
	}

	ctrlCell.appendChild(ctrlButton);
	nameCell.innerHTML = seqName;
	chainCell.innerHTML = chain;
	regionCell.innerHTML = region;
	hiddenSeqCell.innerHTML = seq;
	
	//initialize values
	doc.getElementById("sName").value="";
	doc.getElementById("sVal").value="";
	hiddenSeqCell.style.display="none";
	doc.getElementById("chainGroupSH").checked = true;
	doc.getElementById("regionGroupSV").checked = true;
}

function createDBs()
{
	var dbFileList = opener.document.getElementById("dbFileList").innerHTML;
	var seqList = opener.document.getElementById("seqList").innerHTML;

	var seqTable = document.getElementById("currentSeqs");

	var i;
	if(seqTable != null && seqTable.rows.length > 1)
	{
		for(i = 1; i < seqTable.rows.length; ++i)
		{
			var strInfo = "";
			var seqName = seqTable.rows[i].cells[0].childNodes[0].data;
			var chain = seqTable.rows[i].cells[1].childNodes[0].data;
			var region = seqTable.rows[i].cells[2].childNodes[0].data;
			var seq = seqTable.rows[i].cells[5].childNodes[0].data;
			strInfo = seqName + ":" + chain + ":" + region + ":" + seq;
			if(seqList == null || seqList == "")
				seqList = strInfo;
			else
				seqList += ";" + strInfo;
		}
	}

	opener.document.getElementById("seqList").innerHTML = seqList;


	var fileTable = document.getElementById("currentSeqFiles");

	if(fileTable != null && fileTable.rows.length > 1)
	{
		for(i = 1; i < fileTable.rows.length; ++i)
		{
			var strInfo = "";
			var fileName = fileTable.rows[i].cells[0].childNodes[0].data;
			var chain = fileTable.rows[i].cells[1].childNodes[0].data;
			var region = fileTable.rows[i].cells[2].childNodes[0].data;
			
			strInfo = fileName + ":" + chain + ":" + region;
			if(dbFileList == null || dbFileList == "")
				dbFileList = strInfo;
			else
				dbFileList += ";" + strInfo;
		}
	}

	opener.document.getElementById("dbFileList").innerHTML = dbFileList;

	window.close();
}	


function ClearVales(){
		//reset all values
		location.reload();
		
		
		
}

	

	
	function delRow(item1, item2){
		
		
		if( item2 == 1)
		{
			var theTable = document.all.currentSpectrumTable
			//Call script that delete corresponding file
			var out = theTable.rows[item1].cells[0].firstChild.nodeValue;
			var xmlhttp;
			
			if (window.XMLHttpRequest)
			{// code for IE7+, Firefox, Chrome, Opera, Safari
				xmlhttp=new XMLHttpRequest();
			}
			else
			{// code for IE6, IE5
				xmlhttp=new ActiveXObject("Microsoft.XMLHTTP");
			}
			
			//Pass the object to script file
			xmlhttp.open("POST","remove_file.py",true);
			xmlhttp.send(out);
		}
		else
		{
			var theTable = document.all.ptmtable
		}
						
		theTable.deleteRow(item1);
				
		if ( theTable.rows.length ==  1)
		{
			theTable.deleteRow(0);
			theTable.setAttribute("border","0");
			theTable.setAttribute("style","italic");
			var newRow = theTable.insertRow(theTable.rows.length)
			// give the row its own ID
			newRow.id = newRow.uniqueID
			newRow.align = "left"
			newRow.valign = "center"
			newRow.height ="10"
			newRow.style.fontFamily ="Candara"
			newRow.style.fontSize ="15"
			
			var newCell
					
			// an inserted row has no cells, so insert the cells
			newCell = newRow.insertCell(0)
			// give this cell its own id
			newCell.id = newCell.uniqueID
			// display the row's id as the cell text
			newCell.innerText = "No Spectrum files have been added"
			newCell.align ="Center"
			
		}
	
}



	

