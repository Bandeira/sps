
function basename (path, suffix) {
    // Returns the filename component of the path  
    // 
    // version: 1109.2015
    // discuss at: http://phpjs.org/functions/basename    // +   original by: Kevin van Zonneveld (http://kevin.vanzonneveld.net)
    // +   improved by: Ash Searle (http://hexmen.com/blog/)
    // +   improved by: Lincoln Ramsay
    // +   improved by: djmix
    // *     example 1: basename('/www/site/home.htm', '.htm');    // *     returns 1: 'home'
    // *     example 2: basename('ecra.php?p=1');
    // *     returns 2: 'ecra.php?p=1'
    var b = path.replace(/^.*[\/\\]/g, '');
     if (typeof(suffix) == 'string' && b.substr(b.length - suffix.length) == suffix) {
        b = b.substr(0, b.length - suffix.length);
    }
 
    return b;
}

var ar_ext = ['mzXML', 'mgf'];        // array with allowed extensions
var floatRegExp = /^[\+\-]?(\d+(\.\d*)?|\.?\d+)$/;

function checkName(el, to, sbm) {
	
	
	// - www.marplo.net
	// get the file name and split it to separe the extension
	var name = el.value;
	var ar_name = name.split('.');

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
	}
	else {
		// delete the file name, disable Submit, Alert message
		el.value = '';
		alert('".'+ ar_name[1]+ '" is not an accepted file type.  Must be .mgf or mmzXML');
	}
		  
}
	
	
function addSpecFileRow(form) {
	
	var fileName = form.filess.value;
	var description = form.desc.value;
	var parentMass = form.pmass.value;
	var ionTol = form.itol.value;
	var protease = form.protease.value;
	var cGroup = form.cysteine.value;
	
	alert(fileName);
	//Validate the values
	if(description == null || description == "")
	{
		alert("Please enter a file description");
		form.desc.focus();
		return;
	}
	
	if(fileName == null || fileName == "")
	{
		alert("Please enter a spectrum file");
		form.filess.focus();
		return;
	}
	if (!floatRegExp.test(parentMass)) {
		alert("Parent mass tolerance must be a real number.");
		form.pmass.focus();
		return;
	}
	var fPMTol = parseFloat(parentMass);
	
	//alert("checked stuffA");
	if(fPMTol > 2.5 || fPMTol < 0)
	{
		alert("Invalid parent mass tolerance.  Must be between 0 and 2.5 Da");
		form.pmass.focus();
		return;
	}
	if (!floatRegExp.test(ionTol)) {
		alert("Ion tolerance must be a real number.");
		form.itol.focus();
		return;
	}
	//alert("checked stuffB");
	var fITol = parseFloat(ionTol);
	if(fITol > 1.0 || fITol < 0)
	{
		alert("Invalid ion mass tolerance.  Must be between 0 and 1 Da");
		form.itol.focus();
		return;
	}
	
	//alert("building table!!");
	
	var theTable = document.getElementById("currentSpectrumTable");
	alert(theTable.rows.length);
	theTable.setAttribute("border","1");
	
	//If the table is empty then it only has 1 lines
	//1. 'No spectrum files have been loaded'
	if ( theTable.rows.length == 1)
	{
		//delete No spectrum files have been added
		theTable.deleteRow(0)
	
		var hRow = theTable.insertRow(0);
		
		//insert a new header
		var newCell
					
		// an inserted row has no cells, so insert the cells
		newCell = hRow.insertCell(0)
		// give this cell its own id
		newCell.id = newCell.uniqueID
		// display the row's id as the cell text
		newCell.innerText = "File Name"
		newCell.width = "269"
		newCell.align ="left"
		newCell.bgcolor ="#0000FF"
					
		// reuse cell var for second cell insertion
		newCell = hRow.insertCell(1)
		newCell.id = newCell.uniqueID
		newCell.innerText = "PM Tolerance"
		newCell.width = "145"
					
		// reuse cell var for third cell insertion
		newCell = hRow.insertCell(2)
		newCell.id = newCell.uniqueID
		newCell.innerText = "Ion tolerance"
		newCell.width = "145"
					
		// reuse cell var for fourth cell insertion
		newCell = hRow.insertCell(3)
		newCell.id = newCell.uniqueID
		newCell.innerText = "Protease"
		newCell.width = "145"
					
		// reuse cell var for fifth cell insertion
		newCell = hRow.insertCell(4)
		newCell.width = "50"
					
		newCell = hRow.insertCell(5)
		newCell.id = newCell.uniqueID
		newCell.style.display ="none";
		newCell.width = "145"
		
	}
	var hRow = theTable.insertRow(theTable.rows.length);
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
	iTolCell.innertHTML = ionTol;
	protCell.innerHTML = protease;
	
	//initialize values
	document.all.filess.value="";
	document.all.cysteine.value="C,+57";
	document.all.protease.value="Trypsin";
	document.all.pmass.value="2";
	document.all.itol.value="0.5";
		
	
}
	function launchFileTree()
	{
				
				$('#fileTreeDemo_1').fileTree({ root: 'demo/', script: 'jqueryFileTree.php' }, function(file) { 
					alert(file);
				});
				
				$('#fileTreeDemo_2').fileTree({ root: 'demo/', script: 'jqueryFileTree.php', folderEvent: 'click', expandSpeed: 750, collapseSpeed: 750, multiFolder: false }, function(file) { 
					alert(file);
				});
				
				$('#fileTreeDemo_3').fileTree({ root: 'demo/', script: 'jqueryFileTree.php', folderEvent: 'click', expandSpeed: 750, collapseSpeed: 750, expandEasing: 'easeOutBounce', collapseEasing: 'easeOutBounce', loadMessage: 'Un momento...' }, function(file) { 
					alert(file);
				});
				
				$('#fileTreeDemo_4').fileTree({ root: 'demo/', script: 'jqueryFileTree.php', folderEvent: 'dblclick', expandSpeed: 1, collapseSpeed: 1 }, function(file) { 
					alert(file);
				});
				
	
		
	}
	
	function ClearVales()
	{
		//reset all values
		location.reload();
		//delete spectra folder if it exist
		var xmlhttp;
		if (window.XMLHttpRequest)
		{// code for IE7+, Firefox, Chrome, Opera, Safari
			xmlhttp=new XMLHttpRequest();
		}
		else
		{// code for IE6, IE5
			xmlhttp=new ActiveXObject("Microsoft.XMLHTTP");
		}
		
		xmlhttp.open("POST","clear.py",true);
		xmlhttp.send();
	}

	function loadXMLDoc(item1, item2, item3, item4, item5,item6)
	{
			
			var result = validateForm(item5);
			
			if (result == false )
				return false;
				
			if( item5 == 1)
			{
					var formElement = document.getElementById("firstform");

					//refresh firstform
					
					
					var xmlhttp;
					if (window.XMLHttpRequest)
					{// code for IE7+, Firefox, Chrome, Opera, Safari
						xmlhttp=new XMLHttpRequest();
					}
					else
					{// code for IE6, IE5
						xmlhttp=new ActiveXObject("Microsoft.XMLHTTP");
					}
					
					var Table = document.all.currentSpectrumTable
					var files = document.getElementById('files').files;
					
					//check whether file name is already exist or not
					var i=0; 
					for(i=1; i<Table.rows.length; i++)
					{
						var out = Table.rows[i].cells[0].firstChild.nodeValue;
						if( files[0].name == out)
						{
							alert("Already exist");
							return false;
						}
					}
					
					var i=0; 
					for(i=1; i<Table.rows.length; i++)
					{
						var out = Table.rows[i].cells[5].firstChild.nodeValue;
						if( item6 != out)
						{
							alert("One Cysteine-protecting group must be used");
							return false;
						}
					}
					
					xmlhttp.open("POST","save_file.py",false);
					xmlhttp.send(new FormData(formElement));
					
					addRow(item1, item2, item3, item4, item5, item6);
					
			}
			else
				addRow(item1, item2, item3, item4, item5, item6);
			
	}

	function checkAlphabet(x)
	{
		var chk = true;
		var counter =0;
		var i;
		
		for(i=0; i<= x.length-1 ; i++)
		{
			var ch = x.substring(i,i+1);
			if( (ch>="a" && ch<="z") || (ch>="A" && ch<="Z") ) {
				chk = true;
			}
			else{
				chk=false;
				if (ch=="*") 
				{
					counter = counter+1;
					chk = true;
				}
			}
		}
		
		if( counter > 1)
			chk = false;
			
		if( chk == false )
		{
			alert("Input should be letter,string or *");
			return chk;
		}
	}

	function checkNumber(x)
	{
		var chk = true;
		var counter =0;
		var i;
		
		for(i=0; i<= x.length-1 ; i++)
		{
			var ch = x.substring(i,i+1);
			if( ch >="0" && ch <= "9") {
				chk = true;
			}
			else{
				chk=false;
				if (ch==".") 
				{
					counter = counter+1;
					chk = true;
				}
			}
		}
		
		if( counter > 1)
			chk = false;
		
		if( chk == false )
		{
			alert("Input should be number in range[0, 2.5]");
			return chk;
		}
		
	}
	function validateForm(item5)
	{	
		var chk = false;
		if(item5 == 1)
		{
			var x = document.all.pmass.value;
			var y = document.all.itol.value;

			//Check file updated
			var files = document.getElementById('files').files;
			
			if( files[0] == null)
			{
				alert("Please select file");
				return false;
			}
			
			
		}
		else 
		{
			var x = document.all.mass.value;
			var y = document.all.residue.value;
		}
		
		
		if (x==null || x=="" || y==null || y=="")
		{
			  alert("You must fill out necessary information");
			  return false;
		}
		
		// check illegal input
		if ( checkNumber(x) == false )
			return false;
			
		if( item5 == 1)
		{
			if ( x < 0 || x > 2.5)
			{
				alert("x Input should be number in range[0,2.5]");
				return false;
			}

			if ( checkNumber(y) == false )
				return false;
			
			if ( y < 0 || y > 2.5)
			{
				alert("y Input should be number in range[0,2.5]");
				return false;
			}
		}
		else
		{
			// Check y whether this is alphabet of not
			if( checkAlphabet(y) == false)
				return false;
			
			var selection = document.getElementsByName('rtype');				
		
			for (i=0; i<selection.length; i++)
			{
				if (selection[i].checked == true)
				{
					chk = true;
				}
			}
			
			if( chk == false)
				alert("Please select Type");
				
			return chk;
		}
		

		
		return true;
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

	function addRow(item1, item2, item3, item4, item5, item6) 
	{
		
		if( item5 == 2)
		{
			//validity check
			var result = validateForm(item5);
			
			if( result == false)
				return false;
		
		}
		
		
		if( item2 )
		{
			// assign long reference to shorter var name
			if( item5 == 1 )
			{
				var theTable = document.all.currentSpectrumTable
				var x=document.getElementById('showResult').rows[0].cells;
				//x[1].innerHTML=parseInt(x[1].innerHTML)+1;
				x[1].innerHTML=parseInt(theTable.rows.length);
			}
			else
				var theTable = document.all.ptmtable
						
				theTable.setAttribute("border","1");
		
				if ( theTable.rows.length == 1 )
				{
					//delete No spectrum files have been added
					theTable.deleteRow(0)
					
					//generate new header
					// append new row to the end of the table
					var newRow = theTable.insertRow(theTable.rows.length)
					// give the row its own ID
					newRow.id = newRow.uniqueID
					newRow.align = "left"
					newRow.valign = "center"
					newRow.height ="10"
					newRow.style.fontFamily ="Candara"
					newRow.style.fontSize ="15"
					newRow.setAttribute("bgcolor","#0000FF");
					newRow.style.color= "white";
						
					// declare cell variable
					var newCell
					
					// an inserted row has no cells, so insert the cells
					newCell = newRow.insertCell(0)
					// give this cell its own id
					newCell.id = newCell.uniqueID
					// display the row's id as the cell text
					newCell.innerText = "File Name"
					newCell.width = "269"
					newCell.align ="left"
					newCell.bgcolor ="#0000FF"
					
					// reuse cell var for second cell insertion
					newCell = newRow.insertCell(1)
					newCell.id = newCell.uniqueID
					newCell.innerText = "PM Tolerance"
					newCell.width = "145"
					
					// reuse cell var for third cell insertion
					newCell = newRow.insertCell(2)
					newCell.id = newCell.uniqueID
					newCell.innerText = "Ion tolerance"
					newCell.width = "145"
					
					// reuse cell var for fourth cell insertion
					newCell = newRow.insertCell(3)
					newCell.id = newCell.uniqueID
					newCell.innerText = "Protease"
					newCell.width = "145"
					
					// reuse cell var for fifth cell insertion
					newCell = newRow.insertCell(4)
					newCell.width = "50"
					
					newCell = newRow.insertCell(5)
					newCell.id = newCell.uniqueID
					newCell.style.display ="none";
					newCell.width = "145"
					
					
					
				}
		
				if( item5 == 1 )
					var newRow = theTable.insertRow(theTable.rows.length)
				else
					var newRow = theTable.insertRow(theTable.rows.length-1)
				
				// give the row its own ID
				newRow.id = newRow.uniqueID
				newRow.align = "center"
				newRow.valign = "center"
				newRow.height ="10"
				newRow.style.fontFamily ="Candara"
				newRow.style.fontSize ="15"
				newRow.style.color ="#7F7F7F"
					
				// declare cell variable
				var newCell
				
				// an inserted row has no cells, so insert the cells
				newCell = newRow.insertCell(0)
				// give this cell its own id
				newCell.id = newCell.uniqueID
				// display the row's id as the cell text
				
				// Insert file name 
				if( item5 == 1)
				{
					var files = document.getElementById('files').files;
					newCell.innerText = files[0].name;
				}
				
				newCell.width = "269"
				newCell.align ="left"
				
				
				// reuse cell var for second cell insertion
				newCell = newRow.insertCell(1)
				newCell.id = newCell.uniqueID
				newCell.innerText = item2
				newCell.width = "145"

				
				// reuse cell var for third cell insertion
				newCell = newRow.insertCell(2)
				newCell.id = newCell.uniqueID
				newCell.innerText = item3
				newCell.width = "145"

				
				// reuse cell var for fourth cell insertion
				newCell = newRow.insertCell(3)
				newCell.id = newCell.uniqueID
				
				
				if( item5 == 1)
					newCell.innerText = item4
				else
				{
					var selection = document.getElementsByName('rtype');				
					for (i=0; i<selection.length; i++)
					{
						if (selection[i].checked == true)
						{
							newCell.innerText = selection[i].value;
						}
					}
				}
				newCell.width = "145"
				newCell.align ="left"
				
				// reuse cell var for fifth cell insertion
				newCell = newRow.insertCell(4)
				newCell.width = "50"
				
				var img = document.createElement("img");
				img.src = "images/recycle.jpg";
				img.id = "trash"
				//var ttt = theTable.rows[this.parentNode.parentNode.rowIndex].cell[0];
				img.onclick = function(){delRow(this.parentNode.parentNode.rowIndex, item5)};
				
				newCell.appendChild(img);
				
				
				newCell = newRow.insertCell(5)
				newCell.id = newCell.uniqueID
				newCell.innerText = item6
				newCell.style.display ="none";
				newCell.width = "145"
				
				//initialize values
				document.all.filess.value="";
				document.all.cysteine.value="C,+57";
				document.all.protease.value="Trypsin";
				document.all.pmass.value="2";
				document.all.itol.value="0.5";
				
		}
		
	}


	function loadXMLDocFinal(item1)
	{
			//validate whether directory name is specified
			var projectString = document.all.pdir.value;
			
		
			if( projectString =="")
			{
				alert("Please type project directory path");
				return false;
			}
			
			if( projectString[0]==" ")
			{
				alert("Plese type correct directory path");
				return false;
			}
			
			if( projectString.indexOf("/") != -1)
			{
				alert("Plese type correct directory path");
				return false;
			}
			
			if( projectString.indexOf("'\'") != -1)
			{
				alert("Plese type correct directory path");
				return false;
			}
			
			var xmlhttp;
			if (window.XMLHttpRequest)
			{// code for IE7+, Firefox, Chrome, Opera, Safari
				xmlhttp=new XMLHttpRequest();
			}
			else
			{// code for IE6, IE5
				xmlhttp=new ActiveXObject("Microsoft.XMLHTTP");
			}
		
			xmlhttp.open("POST","folder_exist.py",false);
			xmlhttp.send(item1);
			
		
			xmlhttp.open("POST","fileExist.txt",false);
			xmlhttp.send();
			var folderExist=xmlhttp.responseText;
		
			if( folderExist == "True")
			{
				alert("The folder already exist, please use other name");
				return false;
			}
		
			//return xmlhttp.responseText.indexOf("404 Not Found") > 0 ? false : true;
			
			//check whether the directory is already exist or not
			//if the folder name already exists then return false
						
			var num=0;
			//concatenate files in EXE_DIR/GenoMS base on user check
			var w = document.all.MouseHC.checked;
			var x = document.all.MouseLC.checked;
			var y = document.all.RatHC.checked;
			var z = document.all.RatLC.checked;
			
			
			if( w == true)
				num = num +1;
				
			if( x == true)
				num = num +2;
			
			if( y == true)
				num = num +10;
				
			if( z == true)
				num = num +15;
			
			var theTable = document.all.currentSpectrumTable
			
			
			if(num == 0 || theTable.rows.length==1)
			{
				if( num == 0)
					alert("Missing sequence! Please add/select sequences");
				else
					alert("ERROR: Must specify at least one spectrum file");
				return false;
			}
		
			
			// If "Include common contaminants" checked, copy 
			// commoncontaminants.RS.trie¡¯ and ¡®commoncontaminants.RS.index¡¯ to the project directory
			var common = document.all.cCheck.checked;
			
			if( common == true)
			{
				// Generate directory and copy files in spectra to the directory
				xmlhttp.open("POST","final_fileTwo.py",false);
				xmlhttp.send(item1);
			}
			
			// Generate Mouse_LC,,etc files and "copied files in spectra" to the project directory
			xmlhttp.open("POST","final_file.py",false);
			xmlhttp.send(item1);
			
			//Copy sps.params to tmp.params
			xmlhttp.open("POST","copy_param.py",false);
			xmlhttp.send();	
			
			// Concatenate files based on user input and user sequence, and modify parameter value
			xmlhttp.open("POST","concatenate_file.py",false);
			xmlhttp.send(num);
			
			//Update the project path
			//Update parameter
			xmlhttp.open("POST","param_update.py",false);
			xmlhttp.send("PROJECT_DIR="+item1);
		
			//Insert file names
			var i=0; 
			var max = "=";
			for(i=1; i<theTable.rows.length; i++)
			{
				var out = theTable.rows[i].cells[0].firstChild.nodeValue;
				
				max = max+out+";";
				
			}
			
			//Update parameter
			xmlhttp.open("POST","param_update.py",false);
			xmlhttp.send("INPUT_SPECTRA="+max);
						
			// Find the largest value
			var out = theTable.rows[1].cells[5].firstChild.nodeValue;
				
			//Update parameter
			xmlhttp.open("POST","param_update.py",false);
			xmlhttp.send("FLEXEDMOD="+out);
						
			// Find the correct value
			var i=0; 
			var max1 = 0;
			var max2 = 0;
			var max3 = 0;
			
			for(i=1; i<theTable.rows.length; i++)
			{
				var out = theTable.rows[i].cells[3].firstChild.nodeValue;
				
				if( out == "Trypsin")
					max1 = 1;
				else if( out == "Chymotrypsin")
					max2 = 10;
				else if( out == "Other")
					max3 = 100;	
				else
					max3 = -1;
			}
			var pass;
			if (max1 == 1 && max2 == 0 && max3 == 0)
				pass = "Trypsin";
			else if (max1 == 0 && max2 == 10 && max3 == 0)
				pass = "Chymotrypsin";
			else if(max1 == 0 && max2 == 0 && max3 == 100)
				pass = "Other";
			else 
				pass = "Error"
				
			//Update parameter
			if ( pass != "Error")
			{
				xmlhttp.open("POST","param_update.py",false);
				xmlhttp.send("DIGEST="+pass);
			}
			
			// Find the largest value
			var i=0; 
			var max = 0;
			for(i=1; i<theTable.rows.length; i++)
			{
				var out = theTable.rows[i].cells[1].firstChild.nodeValue;
				
				if( out > max)
				{
					max = out;
				}
			}
			
			//Update parameter
			xmlhttp.open("POST","param_update.py",false);
			xmlhttp.send("TOLERANCE_PM="+max);
			
			
			// Find the largest value
			var i=0; 
			var max = 0;
			for(i=1; i<theTable.rows.length; i++)
			{
				var out = theTable.rows[i].cells[2].firstChild.nodeValue;
				
				if( out > max)
				{
					max = out;
				}
			}
			
			//Update parameter
			xmlhttp.open("POST","param_update.py",false);
			xmlhttp.send("TOLERANCE_PEAK="+max);
			
			//Update parameter if common is on
			if( common == "on")
			{
				xmlhttp.open("POST","param_update.py",false);
				xmlhttp.send("CONTAMINANTS=commoncontaminants.RS.trie");
			}
			else
			{
				xmlhttp.open("POST","param_update.py",false);
				xmlhttp.send("CONTAMINANTS= ");
			}
						
			//Move parameter file to conFile folder
			xmlhttp.open("POST","move_param.py",false);
			xmlhttp.send();
			
			//Copy concFiles to generated directory
			xmlhttp.open("POST","copy_folder.py",true);
			xmlhttp.send(item1);
			
/*			
			//Run executable file
		//	xmlhttp.open("POST","exe_file.py",true);
		//	xmlhttp.send();
*/	
			location.reload();
	}

