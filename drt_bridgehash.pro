; docformat = 'rst'
;
;+
;
; Pass hash to IDL child process.
;
; :Categories:
;    DRTOOLS
;
; :Returns:
;    None
;
; :Params:
;    hash: in, required, type=hash
;      Hash to pass to child process
;    obj: in, required, type=object
;      Object referencing child process
;
; :Keywords:
;    level: in, optional, type=integer, def=-1
;      Level from which to grab hash name to pass; 0 = this routine, 
;      -1 = calling routine, etc.
;    outvarname: in, optional, type=string
;      If set, this is the name for the hash in the child process; otherwise
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
;
; :Copyright:
;    Copyright (C) 2016 David S. N. Rupke
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
pro drt_bridgehash,hashin,obj,level=level,outvarname=outvarname

;  name of temporary variable in child process
   tmpvarname='tmpvar'+string(randomu(!NULL,/ulong),for='(I0)')

;  Get name of hash
   if keyword_set(outvarname) then hashobj=outvarname $
   else begin
      if not keyword_set(level) then level=-1
      hashobj=scope_varname(hashin,level=level)
   endelse
   
   obj->execute, hashobj+"=hash()
   foreach key,hashin.keys() do begin
      obj->setvar, tmpvarname, hashin[key]
      obj->execute, hashobj+'["'+key+'"]='+tmpvarname
   endforeach

   obj->execute, 'undefine, '+tmpvarname

end
