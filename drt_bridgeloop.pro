; docformat = 'rst'
;
;+
;
; Divide execution of a loop among different child processes. Presently does
; not parse output of a child process.
; 
; Structures and hashes are passed to the child process via DRT_BRIDGESTRUCT
; and DRT_BRIDGEHASH.
;
; :Categories:
;    DRTOOLS
;
; :Returns:
;    None.
;
; :Params:
;    nloop: in, required, type=int
;      Number of loop elements.
;    ncores: in, required, type=int
;      Number of cores to split processing over.
;    batchfilein: in, required, type=string
;      Name of text file (can be a normal .PRO) which is parsed and turned into
;      a batch file.
;    batchdirout: in, required, type=string
;      Location where temporary batch files are created to feed to child
;      processes.
;
; :Keywords:
;    invar: in, optional, type=strarr
;      Variables that are to be passed to the child process
;    loopvar: in, optional, type=string, def='i'
;      Name of variable that BATCHFILEIN uses to index the loop.
;    logfiles: in, optional, type=strarr
;      Name of logfile for each child process. If not set, every child process
;      outputs to STDOUT.
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
;      2016sep16, DSNR, created; inspired by R. da Silva's SPLIT_FOR
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
pro drt_bridgeloop,nloop,ncores,batchfilein,batchdirout,$
                   invar=invar,loopvar=loopvar,logfiles=logfiles

   secstr = string(systime(1),format='(I0)')
   if ~keyword_set(loopvar) then loopvar='i'
   if ~keyword_set(logfiles) then logfiles=strarr(ncores)

;  Get number of loop elements in each spawned process
   min_per_proc=floor(double(nloop)/double(ncores))
   loopbeg=indgen(ncores) ; beginning loop index to process on each core
   loopend=loopbeg+(min_per_proc-1d)*ncores ; end loop index on each core
   remainder=round((double(nloop)/double(ncores)-double(min_per_proc))*double(ncores))
   for i=0,remainder-1 do loopend[i]+=ncores ; sprinkle remainder among cores
   loopbeg = string(loopbeg,format='(I0)')
   loopend = string(loopend,format='(I0)')

;  Process input procedure
   openr,batlunin,batchfilein,/get_lun
   outarr=!NULL
   line=''
   while ~ EOF(batlunin) do begin
      readf,batlunin,line
      line = strtrim(line,2) ; get rid of trailing / leading blanks
      iscomment = strmid(line,0,1) ; check if it's a comment
      iscont = strmid(line,0,1,/rev) ; check if there's a continuation
      if iscont eq '$' then tail = '' else tail = ' &$' ; if there's no cont.
      if iscomment ne ';' AND iscomment ne '' then begin
         isproorend = strmid(line,0,3) ; make sure not a pro declaration or end
         length=strlen(line)
         if isproorend ne 'pro' AND (isproorend ne 'end' OR length gt 3) then $
            outarr=[outarr,line+tail]
      endif
   endwhile

;  Write batch files and start up child processes
   obridge=!NULL
   for i=0, ncores-1 do begin

      batchfileout='drt_bridgeloop_'+secstr+'_'+string(i,for='(I0)')
      openw,batlunout,batchdirout+batchfileout+'.pro',/get_lun
      printf,batlunout,'for '+loopvar+'='+loopbeg[i]+','+loopend[i]+','+$
             string(ncores,format='(I0)')+' do begin &$'
      for j=0,n_elements(outarr)-1 do printf,batlunout,outarr[j]
      printf,batlunout,'endfor'
      free_lun,batlunin,batlunout

      obridge=[obridge, obj_new("IDL_IDLBridge",output=logfiles[i])]
;     Check what type of variables are being passed. TMPVAR holds the variable
;     in this procedure; invar[j] is just the name
      if keyword_set(invar) then begin
         for j=0,n_elements(invar)-1 do begin
            tmpvar=scope_varfetch(invar[j],level=-1)
            if isa(tmpvar,'STRUCT') then $
               drt_bridgestruct,tmpvar,obridge[i],outvar=invar[j] $
            else if isa(tmpvar,'HASH') then $
               drt_bridgehash,tmpvar,obridge[i],outvar=invar[j] $
            else $
               obridge[i]->setvar,invar[j],tmpvar
         endfor
      endif
      obridge[i]->execute,'@'+batchdirout+batchfileout,/nowait
   endfor

;  Ping child processes to check completeness. If all cores are done executing,
;  exit
   coresdone=0
   while coresdone lt ncores do begin
      coresdone=0
;     check that core is not executing
      for i=0,ncores-1 do coresdone+=(obridge[i]->status() ne 1)   
   endwhile

;  Delete batch files and child processes
   for i=0, ncores-1 do begin
      batchfileout='drt_bridgeloop_'+secstr+'_'+string(i,for='(I0)')
      file_delete,batchdirout+batchfileout+'.pro',/allow,/quiet
      obj_destroy, obridge[i]
   endfor
   
end
