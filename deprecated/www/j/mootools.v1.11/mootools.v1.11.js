function $merge(){
var mix={};
for(var i=0; i < arguments.length; i++){
for(var property in arguments[i]){
var ap=arguments[i][property];
var mp=mix[property];
if(mp && $type(ap)=='object' && $type(mp)=='object')mix[property]=$merge(mp, ap);
else mix[property]=ap;
}
}
return mix;
};
var $extend=function(){
var args=arguments;
if(!args[1])args=[this, args[0]];
for(var property in args[1])args[0][property]=args[1][property];
return args[0];
};
var $native=function(){
for(var i=0, l=arguments.length; i < l; i++){
arguments[i].extend=function(props){
for(var prop in props){
if(!this.prototype[prop])this.prototype[prop]=props[prop];
if(!this[prop])this[prop]=$native.generic(prop);
}
};
}
};
$native.generic=function(prop){
return function(bind){
return this.prototype[prop].apply(bind, Array.prototype.slice.call(arguments, 1));
};
};
$native(Function, Array, String, Number);
function $chk(obj){
return !!(obj||obj === 0);
};
var Abstract=function(obj){
obj=obj||{};
obj.extend=$extend;
return obj;
};
var Window=new Abstract(window);
var Document=new Abstract(document);
window.xpath=!!(document.evaluate);
if(window.ActiveXObject)window.ie=true;
else if(document.childNodes && !document.all && !navigator.taintEnabled)window.webkit=true;
else if(document.getBoxObjectFor != null)window.gecko=true;
Object.extend=$extend;
if(typeof HTMLElement=='undefined'){
var HTMLElement=function(){};
if(window.webkit)document.createElement("iframe");
HTMLElement.prototype=(window.webkit)? window["[[DOMElement.prototype]]"]:{};
}
HTMLElement.prototype.htmlElement=function(){};
var Class=function(properties){
var klass=function(){
return(arguments[0] !== null && this.initialize && $type(this.initialize)=='function')? this.initialize.apply(this, arguments): this;
};
$extend(klass, this);
klass.prototype=properties;
klass.constructor=Class;
return klass;
};
Class.prototype={
extend: function(properties){
var proto=new this(null);
for(var property in properties){
var pp=proto[property];
proto[property]=Class.Merge(pp, properties[property]);
}
return new Class(proto);
}
};
Class.Merge=function(previous, current){
if(previous && previous != current){
var type=$type(current);
if(type != $type(previous))return current;
switch(type){
case 'function':
var merged=function(){
this.parent=arguments.callee.parent;
return current.apply(this, arguments);
};
merged.parent=previous;
return merged;
case 'object': return $merge(previous, current);
}
}
return current;
};
Array.extend({
forEach: function(fn, bind){
for(var i=0, j=this.length; i < j; i++)fn.call(bind, this[i], i, this);
},
map: function(fn, bind){
var results=[];
for(var i=0, j=this.length; i < j; i++)results[i]=fn.call(bind, this[i], i, this);
return results;
},
contains: function(item, from){
return this.indexOf(item, from)!= -1;
},
extend: function(array){
for(var i=0, j=array.length; i < j; i++)this.push(array[i]);
return this;
}
});
Array.prototype.each=Array.prototype.forEach;
Array.each=Array.forEach;
function $each(iterable, fn, bind){
if(iterable && typeof iterable.length=='number' && $type(iterable)!= 'object'){
Array.forEach(iterable, fn, bind);
} else {
for(var name in iterable)fn.call(bind||iterable, iterable[name], name);
}
};
Array.prototype.test=Array.prototype.contains;
String.extend({
test: function(regex, params){
return(($type(regex)=='string')? new RegExp(regex, params): regex).test(this);
}
});
var Element=new Class({
initialize: function(el, props){
if($type(el)=='string'){
if(window.ie && props &&(props.name||props.type)){
var name=(props.name)? ' name="'+props.name+'"':'';
var type=(props.type)? ' type="'+props.type+'"':'';
delete props.name;
delete props.type;
el='<'+el+name+type+'>';
}
el=document.createElement(el);
}
el=$(el);
return(!props||!el)? el:el.set(props);
}
});
function $(el){
if(!el)return null;
if(el.htmlElement)return Garbage.collect(el);
if(kfm_inArray(el,[window, document]))return el;
var type=$type(el);
if(type=='string'){
el=document.getElementById(el);
type=(el)? 'element':false;
}
if(type != 'element')return null;
if(el.htmlElement)return Garbage.collect(el);
if(kfm_inArray(el.tagName.toLowerCase(),['object', 'embed']))return el;
$extend(el, Element.prototype);
el.htmlElement=function(){};
return Garbage.collect(el);
};
Element.extend=function(properties){
for(var property in properties){
HTMLElement.prototype[property]=properties[property];
Element.prototype[property]=properties[property];
Element[property]=$native.generic(property);
var elementsProperty=(Array.prototype[property])? property+'Elements':property;
}
};
Element.extend({
set: function(props){
for(var prop in props){
var val=props[prop];
this.setProperty(prop, val);
}
return this;
},
setProperty: function(property, value){
var index=Element.Properties[property];
if(index)this[index]=value;
else this.setAttribute(property, value);
return this;
}
});
Element.Properties=new Abstract({
'class': 'className', 'for': 'htmlFor', 'colspan': 'colSpan', 'rowspan': 'rowSpan',
'accesskey': 'accessKey', 'tabindex': 'tabIndex', 'maxlength': 'maxLength',
'readonly': 'readOnly', 'frameborder': 'frameBorder', 'value': 'value',
'disabled': 'disabled', 'checked': 'checked', 'multiple': 'multiple', 'selected': 'selected'
});
var Garbage={
elements: [],
collect: function(el){
if(!el.$tmp){
Garbage.elements.push(el);
el.$tmp={'opacity': 1};
}
return el;
},
trash: function(elements){
for(var i=0, j=elements.length, el; i < j; i++){
if(!(el=elements[i])|| !el.$tmp)continue;
for(var p in el.$tmp)el.$tmp[p]=null;
for(var d in Element.prototype)el[d]=null;
Garbage.elements[Garbage.elements.indexOf(el)]=null;
el.htmlElement=el.$tmp=el=null;
}
}
};
var Event=new Class({
initialize: function(event){
if(event && event.$extended)return event;
this.$extended=true;
event=event||window.event;
this.event=event;
this.type=event.type;
this.target=event.target||event.srcElement;
if(this.target.nodeType==3)this.target=this.target.parentNode;
this.shift=event.shiftKey;
this.control=event.ctrlKey;
this.alt=event.altKey;
this.meta=event.metaKey;
if(this.type.indexOf('key')>-1){
this.code=event.which||event.keyCode;
for(var name in Event.keys){
if(Event.keys[name]==this.code){
this.key=name;
break;
}
}
if(this.type=='keydown'){
var fKey=this.code - 111;
if(fKey > 0 && fKey < 13)this.key='f'+fKey;
}
this.key=this.key||String.fromCharCode(this.code).toLowerCase();
} else if(this.type.test(/(click|mouse|menu)/)){
this.page={
'x': event.pageX||event.clientX+document.documentElement.scrollLeft,
'y': event.pageY||event.clientY+document.documentElement.scrollTop
};
this.client={
'x': event.pageX ? event.pageX - window.pageXOffset:event.clientX,
'y': event.pageY ? event.pageY - window.pageYOffset:event.clientY
};
this.rightClick=(event.which==3)||(event.button==2);
}
return this;
},
stop: function(){
return this.stopPropagation().preventDefault();
},
stopPropagation: function(){
if(this.event.stopPropagation)this.event.stopPropagation();
else this.event.cancelBubble=true;
return this;
},
preventDefault: function(){
if(this.event.preventDefault)this.event.preventDefault();
else this.event.returnValue=false;
return this;
}
});
