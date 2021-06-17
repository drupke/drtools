; docformat = 'rst'
;
;+
;
; Function to compute B-band luminosity from m_B.
; Solar values from Binney and Merrifield Table 2.1.
;
; :Categories:
;    none
;
; :Returns:
;    log( L_B / L_sun )
;
; :Params:
;    appmb: in, required, type=double
;      Apparent B-band magnitude.
;    dl: in, required, type=double
;      Luminosity distance.
;
; :Keywords:
;    absmb: out, optional, type=double
;      Absolute B-band magnitude
; 
; :Author:
;    David S. N. Rupke::
;      Rhodes College
;      Department of Physics
;      2000 N. Parkway
;      Memphis, TN 38112
;      drupke@gmail.com
;
; :History:
;    ChangeHistory::
;      2006oct11, DSNR, created
;      2007may17, DSNR, ported to stand-alone file
;      2015nov05, DSNR, added documentation and copyright
;      2021jun17, DSNR, moved to DRTOOLS
;    
; :Copyright:
;    Copyright (C) 2015--2021 David S. N. Rupke
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
function drt_lb,appmb,dl,absmb=absmb

   ; Compute distance modulus
   distmod=25d + 5d*alog10(dl)
   ; Convert apparent to absolute magnitude
   absmb=appmb-distmod
   ; Compute L_B / L_B_sun using reference to solar M_B
   lb=-0.4d*(absmb-5.48d)
   ; Convert from L_B / L_B_sun to L_B / L_sun using L_B_sun/L_sun
   lb=lb+alog10(4.67d32/3.826d33)

   return,lb
end
