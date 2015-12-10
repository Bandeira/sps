window.kfm_textfile_keybinding=function(e){
	e=new Event(e);
	if(e.code!=27)return;
	e.stopPropagation();
	kfm_textfile_close();
}
