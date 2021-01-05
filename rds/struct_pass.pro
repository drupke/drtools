;+
;PURPOSE
;	to pass structures to idl bridge objects
;SYNTAX
;	struct_pass, str, obj
;INPUTS
;	str: the structure you want to pass
;	obj: the object of bridge you want to pass to
;	level: set to how many levels back you want to get the
;		structure name from... e.g. set to -1 if you want
;		the name from the calling procedure adn -2 if one further back.
;		same syntax as scope_varname
;EXAMPLE	
;  obridge=obj_new("IDL_IDLBridge", output='')
;  a={tag1:0., tag2:findgen(10), tag3:'hello'}
;  a=[a,a]
;
;  struct_pass, a, obridge
;
;  obj_destroy, obridge
;
;DEPENDENCIES
;	struct_addtags
;Written by J. Arnold & R. da Silva, UCSC, 8-31-2010
;-
pro struct_pass, str, obj, level=level
   if not keyword_set(level) then level=-1
   strname=scope_varname(str, level=level)
   fluff='x1s234'
   tagnames=tag_names(str)
   ntags=n_elements(tagnames)
   obj->setvar, 'tg_names_', tagnames
   strname1=fluff+strname ; name of temporary structure
   for j=0,n_elements(str)-1 do begin
;     Create variable named fluff containing jth element of first tag in 
;     structure
      obj->setvar, fluff, str[j].(0)
;     Create structure with the first tag name populated by the jth element
;     of this tag.
      obj->execute, strname1+"=create_struct(tg_names_[0],"+fluff+" )"
;     Do the same for the other tag names
      if ntags GT 1 then begin
         for i=1, ntags-1 do begin
            obj->setvar, fluff, str[j].(i)
;           Create an intermediate structure with the appropriate tag name ...
            obj->execute, fluff+strname1+"=create_struct(tg_names_["+$
                          string(i)+"],"+fluff+" )"
;           ... and then add it to the whole structure.
            obj->execute, strname1+"=struct_addtags("+$
                          strname1+','+fluff+strname1+')'
         endfor
      endif
;     Concatenate the jth element with the previous elements.
      if j EQ 0 then obj->execute, strname+'='+strname1 else $
      obj->execute, strname+'=['+strname+','+strname1+']'
   endfor
   
   if ntags GT 1 then  obj->execute, 'undefine,'+fluff+strname
   obj->execute, 'undefine, '+strname1 
   obj->execute, 'undefine, tg_names_'
   obj->execute, 'undefine, '+fluff
end
