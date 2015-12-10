<?php
require_once('initialise.php');
if(!isset($_POST['name']) || !isset($_POST['value'])) die ('error("post value missing");');
$sn=$_POST['name'];
$uid=$kfm->user_id;
$value=$_POST['value'];
if($kfm->sdef[$sn]['type']=='select_list'){
	if(!isset($_POST['checked']))die ('error("property checked must be given for a select list");');
	$ch=$_POST['checked'];
	$sval=implode(',',$kfm->setting($sn));
	if(isset($_POST['clean']) && $_POST['clean'])$sval='';
	if($ch){
		$sval.=','.$value;
	}else{
		$sval=preg_replace('/'.$value.',*/','',$sval);
	}
	$sval=trim($sval, ' ,');
	$value=$sval;
}
$s=db_fetch_row('SELECT id, usersetting FROM '.KFM_DB_PREFIX.'settings WHERE name="'.$sn.'" and user_id='.$uid);
if(count($s)){
	if(!$s['usersetting'] && $kfm->user_status !=1) die ('error("No rights to change this setting");');	
	$kfmdb->query('UPDATE '.KFM_DB_PREFIX.'settings SET value="'.mysql_escape_string($value).'" WHERE name="'.$sn.'" AND user_id='.$uid);
}else{
	if(!$kfm->isUserSetting($sn) && $kfm->user_status!=1) die ('error("No rights to change this setting");');	
	$kfmdb->query('INSERT INTO '.KFM_DB_PREFIX.'settings (name, value, user_id, usersetting) VALUES ("'.$sn.'","'.mysql_escape_string($value).'", '.$uid.',1)');
	if($kfm->user_id!=1) echo '$("#todefault_'.$sn.'").fadeIn();';
}
if($sn=='theme')$kfm_session->set('theme',$value);
echo 'style_usersetting("'.$sn.'");';
?>
message('setting changed');
