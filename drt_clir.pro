; docformat = 'rst'
;
;+
;
; Function to compute infrared luminosity:
; L(8-1000um) from Sanders and Mirabel 1996, ARA&A
;
; :Categories:
;    DRTOOLS
;
; :Returns:
;    log[ L_IR / L_sun ]
;
; :Params:
;    f12: in, required, type=double
;       12um IRAS flux density, in Jy
;    f25: in, required, type=double
;       25um IRAS flux density, in Jy
;    f60: in, required, type=double
;       60um IRAS flux density, in Jy
;    f100: in, required, type=double
;       100um IRAS flux density, in Jy
;    dist: in, required, type=double
;       luminosity distance, in Mpc
;
; :Keywords:
;    inlog: in, optional, type=boolean
;       if input flux densities are log(f_nu)
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
;       2006feb16, DSNR, created
;       2006nov22, DSNR, ported to stand-alone file
;       2021jun17, DSNR, moved to DRTOOLS
;
; :Copyright:
;    Copyright (C) 2021 David S. N. Rupke
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
function drt_clir,f12,f25,f60,f100,dist,inlog=inlog

   if keyword_set(inlog) then begin
      f12 = 10d^f12
      f25 = 10d^f25
      f60 = 10d^f60
      f100 = 10d^f100
   endif
   ; log(d_L/cm)
   distcm = alog10(dist) + alog10(3.08567d) + 24d
   ; log(f_IR / W m^-2)
   fir = alog10( 1.8d * (13.48d*f12 + 5.16d*f25 + 2.58d*f60 + f100) ) - 14
   ; log(f_IR / erg s^-1 cm^-2):  m^-2 -> cm^-2 is -4d, W --> erg is +7d
   fir = fir + 3d
   ; log(f_IR / L_sun)
   fir = fir - alog10(3.826d) - 33d
   ; log(L_IR / L_sun)
   lir = fir + alog10(4d*!DPI) + 2d*distcm

   return,lir
end
