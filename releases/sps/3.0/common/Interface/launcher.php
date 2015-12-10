<html>

<head>
	
	<title>Valens Job List</title>
	<link href="styles/main.css" rel="stylesheet" type="text/css" />
    <script src="jquery.js" type="text/javascript"></script>	
</head>
	<body>
<div id="bodywrapper">
		<div class="firstform" id="firstForm">
		<table class="firstform" id="tableForm">
			<tr>
			  <th align="left">Current Project Status</th>
			</tr>
			<tr>
			  <td>
			  <table align="center" width="100%" id="status" name="status" class="ptm">
			    <tr>
			      <th>Project Name</th>
			      <th>Launch Time</th>
			      <th>Status</th>
                              <th>Results</th>
			    </tr>
			    <tr>
					<td><?php echo basename($_GET['projectDir']) ?></td>		
					<td id="launchtime"></td>	    
					<td id="projectstatus">Launching</td>
                                        <td id="results"></td>
			    </tr>
			
			  </table>
			  </td>
			</tr>
		</table>
		</div>
		<script type="text/javascript">
			paramFile= '<?php echo $_GET['paramFile'] ?>';
			exeDir =   '<?php echo $_GET['exeDir'] ?>';
			projectDir =  '<?php echo $_GET['projectDir'] ?>';
//docRoot = '<?php $f=$_SERVER["DOCUMENT_ROOT"];echo(realpath($f)); ?>';
//docRoot = '<?php echo $_SERVER["DOCUMENT_ROOT"]; ?>';
docRoot='<?php echo dirname(__FILE__);?>';
		count = 0;


			params = 'paramFile=' + paramFile + '&exeDir=' + exeDir + '&projectDir=' + projectDir;
			statusTxt = "Starting Program";
			processing = false;
                        reportURL = projectDir + "/report/index.html";

index = reportURL.indexOf(docRoot);

reportURL = reportURL.substring(index+docRoot.length+1,reportURL.length);

//checks the status of the running program
			function checkstatus(){
				if (!processing)
       			{
					timeout = setTimeout(function(){
       					 
					processing=true;
					$.ajax({
						url: 'checkstatus.php',
						data: params,
						error: function(jqXHR, textStatus, errorThrown){
							$('#projectstatus').html("JS Error:" + textStatus + "<br />" + errorThrown);
							processing=false;
						},
						success: onComplete
					});
				
					}, 15000);
				}
				return;
			}
			function onComplete(data){
				processing=false;
				successTxt = data;
				$('#projectstatus').html(successTxt);
				
				if(successTxt.indexOf("Finished") != 0 && successTxt.indexOf("Error") != 0)
				  checkstatus();
				else if(successTxt.indexOf("Finished") == 0)
				  $('#results').html("<a href=" + reportURL + ">Reports</a>");
			}
			//This function waits till everything is loaded.
		  $(function() {
				// Handler for .ready() called.
				//Start program
				$('#launchtime').html(Date());
				$('#projectstatus').html(status);
				$.ajax({
					url: 'executeCSPS.php',
					data: params,
					error: function(jqXHR, textStatus, errorThrown){
						$('#projectstatus').html("JS Error:" + textStatus + "<br />" + errorThrown);
						
					},
					success: function(data){
						statusTxt = data;
						$('#projectstatus').html(statusTxt);
						 checkstatus();
					}
				});
				
				
			});
		 
		</script>
</div>
	</body>
</html>

