; docformat = 'rst'
;
;+
;
; Compute electron density using recipe of Sanders, Shapley et al. 2016.
; 
; Density (+/- errors) are set to maximum or minimum values if it exceeds a threshold.
;
; :Categories:
;    DRTOOLS
;    
; :Returns:
;    N-element array, or Nx3 array if errors, of log(density/cm^-2).
;
; :Params:
;    rat: in, requipred, type=dblarr(N)
;      Line ratios
;    type: in, required, type=string
;      's2' for [SII] ratio, 'o2' for [OII]
;      
; :Keywords:
;    denerrlo: out, optional, type=dblarr(N)
;      Output density lower error.
;    denerrhi: out, optional, type=dblarr(N)
;      Output density upper error.
;    err: in, optional, type=err
;      if present, line ratio errors for computing error in density. Assumed to
;      be linear and one-sided in line ratio. Ouptut errors are two-sided in
;      the log, first low and then high.
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
;      2021dec16, DSNR, created
;      2022jan28, DSNR, error now outputs through keyword
;      2022may24, DSNR, changed minrat, maxrat
;
; :Copyright:
;    Copyright (C) 2021-2022 David S. N. Rupke
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
function drt_nebden,rat,type,err=err,denerrlo=denerrlo,denerrhi=denerrhi
   bad = 1d99
   ;  Values for computing electron density
   ;  from Sanders, Shapley, et al. 2016
   ; [OII] and [SII], respectively
   minrat = [0.3839d,0.4375d]
   maxrat = [1.4558d,1.4484d]
   a = [0.3771d,0.4315d]
   b = [2468d,2107d]
   c = [638.4d,627.1d]
   ;maxden = alog10((c*minrat - a*b)/(a - minrat))
   ;minden = alog10((c*maxrat - a*b)/(a - maxrat))
   ; According to Sanders et al., the min/max ratios correspond to these limiting values
   ; We conservatively set the actual limits factor of 10 higher and lower, respectively
   minden = [1.,1.]
   maxden = [4.,4.]
   if type eq 'o2' then i=0 else i=1
   den = alog10((c[i]*rat - a[i]*b[i])/(a[i] - rat))
   ; Errors
   if keyword_set(err) then begin
      ; from delt_rat = drho/drat * delt_rat
      ;denerr = abs((c[i]/(c[i]*rat - a[i]*b[i]) + 1d/(a[i] - rat))/alog(10d)) * err
      ; from simply adding or subtracting the ratio
      denerrlo = alog10((c[i]*(rat+err) - a[i]*b[i])/(a[i] -(rat+err)))
      denerrhi = alog10((c[i]*(rat-err) - a[i]*b[i])/(a[i] -(rat-err)))
      ;  Check for hitting upper/lower ratio limit. Set 
      ilo = where(rat+err ge maxrat[i] OR denerrlo le minden[i],ctlo)
      ihi = where(rat-err le minrat[i] OR (denerrhi ge maxden[i] AND denerrhi ne bad),cthi)
      if ctlo gt 0 then denerrlo[ilo] = minden[i]
      if cthi gt 0 then denerrhi[ihi] = maxden[i]
   endif
   ; Error bar is difference with measured value. Difference for limit is 0
   ; according to formula above.
   denerrlo = den - denerrlo
   denerrhi = denerrhi - den
;  Check for hitting upper/lower ratio limit
   ilo = where(rat ge maxrat[i] OR den le minden[i],ctlo)
   ihi = where(rat le minrat[i] OR (den ge maxden[i] AND den ne bad),cthi)
   if ctlo gt 0 then begin
      den[ilo] = minden[i]
      denerrlo[ilo] = 0d
   endif
   if cthi gt 0 then begin
      den[ihi] = maxden[i]
      denerrhi[ihi] = 0d
   endif

   return,den
end