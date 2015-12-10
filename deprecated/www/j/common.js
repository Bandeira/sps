// see ../license.txt for licensing
function _(str,context,vars,htmlonly){ // translations
	var el;
	if(typeof(str)=='array' || typeof(str)=='object'){
		context=str[1];
		vars=str[2];
		str=str[0];
	}
	if(typeof(str)!='string')alert(typeof(str)+"\n"+str);
	// { set up variables
	if(!context)context='kfm';
	if(!vars)vars=[];
	if(!window.lang)window.lang=[];
	if(!window.lang[kfm_vars.lang])window.lang[kfm_vars.lang]=[];
	// }
	if(htmlonly){
		if(window.lang[kfm_vars.lang][str]){
			el='<span class="kfmlang kfmlang_'+str.toLowerCase().replace(/[^a-z0-9]/g,'')+'">'+window.lang[kfm_vars.lang][str].text+'</span>';
			window.lang[kfm_vars.lang][str].requested=1;
		}
		else el='<span class="kfmlang kfmlang_'+str.toLowerCase().replace(/[^a-z0-9]/g,'')+'">...</span>';
	}
	else{
		el=document.createElement('span');
		el.className='kfmlang kfmlang_'+str.toLowerCase().replace(/[^a-z0-9]/g,'');
		el.original=[str,context,vars];
		if(window.lang[kfm_vars.lang][str]){
			el.appendChild(document.createTextNode(window.lang[kfm_vars.lang][str].text));
			window.lang[kfm_vars.lang][str].requested=1;
		}
		else el.appendChild(document.createTextNode('...'));
	}
	if(!window.lang[kfm_vars.lang][str]){
		x_kfm_translate(str,context,kfm_vars.lang,function(res){
			if(window.lang[kfm_vars.lang][str])return; // in case the string is requested multiple times
			window.lang[kfm_vars.lang][str]={
				'text':res,
				'requested':1
			}
			$j('.kfmlang_'+str.toLowerCase().replace(/[^a-z0-9]/g,'')).text(window.lang[kfm_vars.lang][str].text);
		});
	}
	return el;
}
function kfm_addWidget(obj){
	kfm_widgets.push(obj);
}
function array_remove_values(arr,vals){
	if($type(vals)!='array')vals=[vals];
	var i,tmp=[];
	for(i=0;i<arr.length;++i)if(vals.indexOf(arr[i])==-1)tmp.push(arr[i]);
	return tmp;
}
function clearSelections(){
	window.getSelection().removeAllRanges();
}
function getOffset(el,s) {
	if(!el){
		return 0;
	}
	var n=parseInt(el['offset'+s],10),p=el.offsetParent;
	if(p){
		n+=getOffset(p,s)-parseInt(p['scroll'+s],10);
	}
	return n;
}
function getWindowScrollAt(){
	return {x:window.pageXOffset,y:window.pageYOffset};
}
function kfm_kaejax_do_call(func_name,args){
	var uri=function_urls[func_name];
	if(!window.kfm_kaejax_timeouts[uri]){
		window.kfm_kaejax_timeouts[uri]={t:setTimeout('kfm_kaejax_sendRequests("'+uri+'")',1),c:[],callbacks:[]};
	}
	var l=window.kfm_kaejax_timeouts[uri].c.length,v2=[];
	for(var i=0;i<args.length-1;++i){
		v2[v2.length]=args[i];
	}
	window.kfm_kaejax_timeouts[uri].c[l]={f:func_name,v:v2};
	window.kfm_kaejax_timeouts[uri].callbacks[l]=args[args.length-1];
}
function kfm_kaejax_sendRequests(uri){
	var t=window.kfm_kaejax_timeouts[uri],callbacks=window.kfm_kaejax_timeouts[uri].callbacks;
	t.callbacks=null;
	window.kfm_kaejax_timeouts[uri]=null;
	var x=new XMLHttpRequest(),post_data="kaejax="+escape(Json.toString(t)).replace(kfm_regexps.plus,'%2B').replace(kfm_regexps.ascii_stuff,'%u00$1').replace(/\n/g,' ');
	post_data=kfm_sanitise_ajax(post_data);
	x.open('POST',uri,true);
	x.setRequestHeader("Method","POST "+uri+" HTTP/1.1");
	x.setRequestHeader("Content-Type","application/x-www-form-urlencoded");
	x.onreadystatechange=function(){
		if(x.readyState!=4){
			return;
		}
		var r=x.responseText;
		if(r.substring(0,5)=='error'){
			return kfm.alert(r);
		}
		var v=eval('('+unescape(r)+')');
		var f,p,i;
		if(v.errors.length)kfm.showErrors(v.errors);
		if(v.messages.length)kfm.showMessages(v.messages);
		for(i=0;i<t.c.length;++i){
			f=callbacks[i];
			p=[];
			if($type(f)=='array'){
				p=f;
				f=f[0];
			}
			f(v.results[i],p);
		}
	};
	x.send(post_data);
}
function loadJS(url,id,lang,onload){
	var i=0,el;
	for(;i<loadedScripts.length;++i){
		if(loadedScripts[i]==url){
			return 0;
		}
	}
	loadedScripts.push(url);
	el=document.createElement('script');
	el.type='text/javascript';
	if(id){
		el.id=id;
	}
	if(lang){
		el.lang=lang;
	}
	if(kfm_kaejax_is_loaded&&/\.php/.test(url)){
		url+=(/\?/.test(url)?'&':'?')+'kfm_kaejax_is_loaded';
	}
	if(onload){
		el.onload_triggered=0;
		el.onload=function(){
			if(!this.onload_triggered++){
				eval(onload);
			}
		};
		el.onreadystatechange=function(){
			if(this.readyState=='loaded'||this.readyState=='complete'){
				if(!this.onload_triggered++){
					eval(onload);
				}
			}
		};
	}
	el.src=url;
	document.getElementsByTagName('HEAD')[0].appendChild(el);
	return 1;
}
function newInput(n,t,v,cl){
	var b;
	if(!t){
		t='text';
	}
	switch(t){
		case 'checkbox': // {{{
			b=document.createElement('input');
			b.type='checkbox';
			b.style.width='auto';
			break;
		// }}}
		case 'textarea': // {{{
			b=document.createElement('textarea');
			break;
		// }}}
		default: // {{{
			b=document.createElement('input');
			b.type=t;
		// }}}
	}
	b.id=n;
	b.name=n;
	if(v){
		if(t=='checkbox'){
			$extend(b,{checked:'checked',defaultChecked:'checked'});
		}
		else b.value=v;
	}
	if(cl)b.className=cl;
	return b;
}
function newLink(h,t,id,c,title){
	if(!title)title='';
	var a=document.createElement('a');
	a.id=id;
	a.className=c;
	a.href=h;
	a.title=title;
	a.innerHTML=t;
	return a;
}
function newSelectbox(name,keys,vals,s,f){
	var el2,el3,s2=0,i=0;
	el2=document.createElement('select');
	el2.id=name;
	if(!s){
		s=0;
	}
	if(!vals){
		vals=keys;
	}
	for(;i<vals.length;++i){
		var v1=vals[i].toString();
		var v2=v1.length>20?v1.substr(0,27)+'...':v1;
		el3=document.createElement('option');
		el3.value=keys[i];
		el3.title=v1;
		el3.innerHTML=v2;
		if(keys[i]==s)s2=i;
		el2.appendChild(el3);
	}
	el2.selectedIndex=s2;
	if(f)el2.onchange=f;
	return el2;
}
function newText(a){
	return document.createTextNode(a);
}
function $defined(obj){
	return (obj != undefined);
}
function $pick(obj, picked){
	return $defined(obj) ? obj : picked;
}
function $type(obj){
	if (!$defined(obj)) return false;
	if (obj.htmlElement) return 'element';
	var type = typeof obj;
	if (type == 'object' && obj.nodeName){
		switch(obj.nodeType){
			case 1: return 'element';
			case 3: return (/\S/).test(obj.nodeValue) ? 'textnode' : 'whitespace';
		}
	}
	if (type == 'object' || type == 'function'){
		switch(obj.constructor){
			case Array: return 'array';
			case RegExp: return 'regexp';
			case Class: return 'class';
		}
		if (typeof obj.length == 'number'){
			if (obj.item) return 'collection';
			if (obj.callee) return 'arguments';
		}
	}
	return type;
}
if(window.ie){
	XMLHttpRequest=function(){
		var l=(ScriptEngineMajorVersion()>=5)?"Msxml2":"Microsoft";
		return new ActiveXObject(l+".XMLHTTP")
	}
	loadJS('j/browser-specific.ie.js');
}
if(window.webkit){
	loadJS('j/browser-specific.konqueror.js');
}
var Json = {
toString: function(obj){
switch($type(obj)){
case 'string':
return '"' + obj.replace(/(["\\])/g, '\\$1') + '"';
case 'array':
return '[' + obj.map(Json.toString).join(',') + ']';
case 'object':
var string = [];
for (var property in obj) string.push(Json.toString(property) + ':' + Json.toString(obj[property]));
return '{' + string.join(',') + '}';
case 'number':
if (isFinite(obj)) break;
case false:
return 'null';
}
return String(obj);
},
evaluate: function(str, secure){
return (($type(str) != 'string') || (secure && !str.test(/^("(\\.|[^"\\\n\r])*?"|[,:{}\[\]0-9.\-+Eaeflnr-u \n\r\t])+?$/))) ? null : eval('(' + str + ')');
}
};
