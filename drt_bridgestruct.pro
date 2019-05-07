; docformat = 'rst'
;
;+
;
; Pass structure to IDL child process.
;
; :Categories:
;    DRTOOLS
;
; :Returns:
;    None
;
; :Params:
;    struct: in, required, type=structure
;      Structure to pass to child process
;    obj: in, required, type=object
;      Object referencing child process
;
; :Keywords:
;    level: in, optional, type=integer, def=-1
;      Level from which to grab structure name to pass; 0 = this routine, 
;      -1 = calling routine, etc.
;    outvarname: in, optional, type=string
;      If set, this is the name for the structure in the child process; otherwise
;      it is pulled from the calling routine.
;
; :Author:
;    David S. N. Rupke::
;      Rhodes College
;      Department of Physics
;      2000 N. Parkway
;      Memphis, TN 38104
;      drupke@gmail.com
;
; :History:
;    ChangeHistory::
;      2016sep16, DSNR, created; inspired by R. da Silva's STRUCT_PASS
;      2019may07, DSNR, changed call to RANDOMU for compatibility with versions
;                       prior to IDL 8.2.2
;
; :Copyright:
;    Copyright (C) 2016--2019 David S. N. Rupke
;
;    This program is free software: you can redistribute it and/or
;    modify it under the terms of the GNU General Public License as
;    published by the Free Software Foundation, either version 3 of
;    the License or any later version.
;
;    This program is distributed in the hope that it will be useful,
;    but WITHOUT ANY WARRANTY; without even the implied warranty of
;    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
;    General Public License for more details.
;
;    You should have received a copy of the GNU General Public License
;    along with this program.  If not, see
;    http://www.gnu.org/licenses/.
;
;-
pro drt_bridgestruct,struct,obj,level=level,outvarname=outvarname

;  name of temporary variable in child process
;   tmpvarname='tmpvar'+string(randomu(!NULL,/ulong),for='(I0)')
   tmpvarname='tmpvar'+string(abs(randomu(!NULL,/long)),for='(I0)')

;  Get name of structure
   if keyword_set(outvarname) then structobj=outvarname $
   else begin
      if not keyword_set(level) then level=-1
      structobj=scope_varname(struct,level=level)
   endelse

   tags=tag_names(struct)
   ntags=n_elements(tags)

;  Loop through elements of array of structures, if there is more than one.   
   for i=0,n_elements(struct)-1 do begin
;     Variable to hold command that will define structure
      structdefstr = '{'
;     Loop through tags.
      for j=0, ntags-1 do begin
;        Variable containing ith element of jth tag in array of structures
         jtmpvarname=tmpvarname+'_'+string(j,format='(I0)')
         structelement=struct[i].(j)
         if isa(structelement,'STRUCT') then begin
            drt_bridgestruct,structelement,obj,outvar=jtmpvarname
         endif else if isa(structelement,'HASH') then begin
            drt_bridgehash,structelement,obj,outvar=jtmpvarname
         endif else begin
            obj->setvar,jtmpvarname,structelement
         endelse
         structdefstr+=tags[j]+':'+jtmpvarname
         if j ne ntags-1 then structdefstr+=','
      endfor
      structdefstr+='}'
;     Execute structure definition.
      if i eq 0 then obj->execute, structobj+'='+structdefstr $
      else obj->execute, structobj = [structobj,structdefstr]
   endfor

;  Clear temporary variables from child process
   for j=0, ntags-1 do begin
      jtmpvarname=jtmpvarname+string(j,format='(I0)')
      obj->execute, 'undefine, '+jtmpvarname
   endfor

end
