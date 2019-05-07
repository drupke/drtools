; docformat = 'rst'
;
;+
;
;  Compute luminosity from flux. Output is in Solar units.
;
; :Categories:
;    SPECTRA
;
; :Returns:
;
; :Params:
;    flux: in, required, type=dblarr
;       Flux in W/m^2. 1 W = 10^7 erg/s, and 1 m^2 = 10^4 cm^2.
;       So 1 W/m^2 = 10^7 erg/s / 10^4 cm^2 = 10^3 erg/s/cm^2
;       If input is in Jy (=10^-26 W/m^2/Hz), then output is 10^26 W/Hz.
;    dist: in, required, type=double
;       Luminosity distance in Mpc.
;
; :Keywords:
;    ergs: in, optional, type=byte
;      Output in erg/s.
;    err: in, optional, type=dblarr
;      Compute error in luminosity; set this keyword to the error array. Output
;      is then an N x 2 array, where N is the number of fluxes, and the second
;      dimension is lum. and error. Doesn't function yet if input is logged.
;    log: in, optional, type=byte
;      Input flux and output luminosity are in log_10 space.
;    z: in, optional, type=double
;      Use redshift instead of distance. D_L then computed using default
;      cosmology in ASTROLIB routine LUMDIST.
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
;      2008jan28, DSNR, created
;      2016sep01, DSNR, added copyright and documentation
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
function drt_linelum,flux,dist,ergs=ergs,err=err,log=log,z=z,watts=watts

   if keyword_set(z) then dist=lumdist(z)

   if keyword_set(log) then begin
      lum = flux + 7d + alog10(4d*!DPI) + 2d*alog10(dist) + 44d + alog10(9.5213d)
      if keyword_set(watts) then lum -= 7d
      if ~ keyword_set(ergs) AND ~ keyword_set(wmsq) $
         then lum = lum - alog10(3.826d) - 33d
   endif else begin
      lum = flux * 1d7 * 4d * !DPI * dist^2d * 9.5213d44
      if keyword_set(watts) then lum /= 1d7
      if keyword_set(err) then begin
         lumerr = err * 1d7 * 4D * !DPI * dist^2d * 9.5213d44
         if keyword_set(watts) then lumerr /= 1d7
      endif
      if ~ keyword_set(ergs) AND ~ keyword_set(watts) then begin
         lum /= 3.826d33
         if keyword_set(err) then lumerr /= 3.826d33
      endif
   endelse

   if ~ keyword_set(err) OR keyword_set(log) then return,lum $
   else return,[[lum],[lumerr]]

end
