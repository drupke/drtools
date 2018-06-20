; docformat = 'rst'
;
;+
;
; Compute the Gaussian function using the parameterization that's in
; MPFITPEAK.
;
; A[0] = peak value
; A[1] = profile center
; A[2] = sigma
; A[3] = constant offset term [optional]
; A[4] = linear offset term [optional]
; 
; :Categories:
;
; :Returns:
;    Scalar or array with Gaussian function evaluated at the input points.
;
; :Params:
;    x: in, required, type=double or dblarr(N)
;      Coordinates at which to evaluate the function.
;    a: in, required, type=dblarr(M)
;      Parameters, as given above.
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
;      2018feb15, DSNR, created
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
function drt_gaussian,x,a

  y = gaussian(x,a,/double)
  if n_elements(a) eq 5 then y = y + a[4]*x

  return,y

end
