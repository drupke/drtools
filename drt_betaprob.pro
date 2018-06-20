; docformat = 'rst'
;
;+
;
; Compute confidence interval from beta distribution. Outputs fractional range
; of confidence interval, normalized to total number.
; Routine from E. Cameron 2011, PASP
;
; :Categories:
;    DRTOOLS
;
; :Returns:
;    Postscript plots.
;
; :Params:
;    k: in, required, type=double
;       Number in which quantity is measured.
;    n: in, required, type=double
;       Total number.
;       
; :Keywords:
;    c: in, optional, type=double, default=0.683
;       Confidence interval for which to compute probability.
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
;      2018jun15, DSNR, created
;
; :Copyright:
;    Copyright (C) 2018 David S. N. Rupke
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
function drt_betaprob,k,n,c=c

   if not keyword_set(c) then c=0.683d

   z = DINDGEN(10000)*0.0001d
   Beta = IBETA(k+1,n-k+1,z)
   ill = VALUE_LOCATE(Beta,(1-c)/2)
   iul = VALUE_LOCATE(Beta,1-(1-c)/2)
   plower = z[ill]
   pupper = z[iul]

   return,[plower,pupper]

end
