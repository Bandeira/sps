window.kfm_files_dragToDirectory=function(e){
	dir_over=e.targetElement.node_id;
	var mv='x_kfm_moveFiles(['+selectedFiles.join(',')+'],'+dir_over+',function(e){if($type(e)=="string")return alert(kfm.lang.CouldNotMoveFiles);kfm_removeFilesFromView(['+selectedFiles.join(',')+'])});kfm_selectNone()';
	var cp='x_kfm_copyFiles(['+selectedFiles.join(',')+'],'+dir_over+',kfm_showMessage);kfm_selectNone()';
	if(kfm_vars.files.drags_move_or_copy == 1)eval(mv); // always Move
	else if(kfm_vars.files.drags_move_or_copy == 2)eval(cp); // always Copy
	else{
		var links=[];
		//links.push([cp,kfm.lang.CopyFiles]);
		context_categories['edit'].add({
			name:'files_copy',
			title:kfm.lang.CopyFiles,
			category:'edit',
			doFunction:function(){eval(cp);}
		});
		//links.push([mv,kfm.lang.MoveFiles,0,!kfm_vars.permissions.file.mv]);
		context_categories['edit'].add({
			name:'files_move',
			title:kfm.lang.MoveFiles,
			category:'edit',
			doFunction:function(){eval(mv)}
		});
		kfm_createContextMenu(e.page,false);
	}
}
