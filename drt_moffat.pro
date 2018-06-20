; docformat = 'rst'
;
;+
;
; Compute the Moffat function using the parameterization that's in
; MPFITPEAK. The model function is
;
; y = A[0]/(u^2 + 1)^A[3] + A[4] + A[5]*x
;
; where
; u = (x - A[1])/A[2]
; A[0] = peak value
; A[1] = profile center
; A[2] = alpha = (FWHM/2)/sqrt[2^(1/beta)-1]
; A[3] = beta = "Moffat index"
; A[4] = constant offset term [optional]
; A[5] = linear offset term [optional]
; 
; :Categories:
;
; :Returns:
;    Scalar or array with Moffat function evaluated at the input points.
;
; :Params:
;    x: in, required, type=double or dblarr(N)
;      Coordinates at which to evaluate the Moffat function.
;    a: in, required, type=dblarr(M)
;      Moffat function parameters, as given above.
;
; :Keywords:
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
;      2013may16, DSNR, created
;      2015may11, DSNR, added documentation, license, and copyright;
;                       corrected documentation to note that
;                       alpha != HWHM exactly.
;    
; :Copyright:
;    Copyright (C) 2015 David S. N. Rupke
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
function drt_moffat,x,a

  u = (x - a[1])/a[2]
  mof = a[0]/(u^2d + 1d)^(a[3])
  if n_elements(a) ge 5 then mof = mof + a[4]
  if n_elements(a) eq 6 then mof = mof + a[5]*x

  return,mof

end
