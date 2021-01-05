;+
;PURPOSE
;  to pass hashes to idl bridge objects
;SYNTAX
;  hash_pass, hash, obj
;INPUTS
;  hash: the hash you want to pass
;  obj: the object of bridge you want to pass to
;  level: set to how many levels back you want to get the
;         hash name from... e.g. set to -1 if you want
;         the name from the calling procedure adn -2 if one further back.
;         same syntax as scope_varname
;EXAMPLE
;
;DEPENDENCIES
;Written by J. Arnold & R. da Silva, UCSC, 8-31-2010
;-
pro hash_pass, hashin, obj, level=level
   if not keyword_set(level) then level=-1
   hashname=scope_varname(hashin, level=level)
   fluff='x1s234'
   keys_list=hashin.keys()
   keys=keys_list.toarray()
   nkeys=hashin.count()
   obj->setvar, 'keys_',keys
   hashname=fluff+hashname
   obj->execute, hashname+"=hash()
   for i=0,nkeys-1 do begin
      obj->setvar, fluff, hashin[keys[i]]
      obj->execute, hashname'='+hashname else $
   endforeach

   if ntags GT 1 then  obj->execute, 'undefine,'+fluff+strname
   obj->execute, 'undefine, '+strname1
   obj->execute, 'undefine, tg_names_'
   obj->execute, 'undefine, '+fluff
end
