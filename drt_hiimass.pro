; docformat = 'rst'
;
;+
;
;  Compute ionized gas mass from Halpha luminosity and electron density.
;
; :Categories:
;    SPECTRA
;
; :Returns:
;
; Mass in solar masses. N-element array if no input errors; Nx2 or Nx3 array 
; otherwise, with elements 2 or 2--3 of second dimension being errors. One-sided
; errors if linear, two-sided otherwise.
;
; :Params:
;    lum: in, required, type=dblarr(N)
;    den: in, required, type=dblarr(N)
;
; :Keywords:
;    errlum: in, optional, type=dblarr(N) or dblarr(N,2)
;      Input luminosity error. One-sided if linear, two-sided if log.
;    errden: in, optional, type=dblarr(N) or dblarr(N,2)
;      Input density error. One-sided if linear, two-sided if log.
;    log: in, optional, type=byte
;      Input luminosity, density, and errors, and output mass, are in log space.
;    masserr: out, optional, type=dblarr(N) or dblarr(N,2)
;      Output error in mass. One-sided if linear, two-sided if log.
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
;      2021dec20, DSNR, created
;      2022jan28, DSNR, error now outputs through keyword
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
function drt_hiimass,lum,den,errlum=errlum,errden=errden,log=log,masserr=masserr

   ; mass per particle, in solar masses
   mumpsm = 1.4d*1.672649d-24 / 1.989d33
   lmumpsm = alog10(1.4d)+alog10(1.672649)-24d -alog10(1.989d)-33d
   ; volume emissivity of Ha = product of recomb. coeff. and photon energy
   ; units erg cm^3 s^-1
   volemis = 2.63d-25
   lvolemis = alog10(2.63d)-25d

   if not keyword_set(log) then begin
      mass = mumpsm * lum / volemis / den
      if keyword_set(errlum) then masserr1 = errlum/lum else masserr1=0d
      if keyword_set(errden) then masserr2 = errden/den else masserr2=0d
      masserr = mass * sqrt(masserr1^2d + masserr2^2d)
      return, mass
   endif else begin
      mass = lmumpsm + lum - lvolemis - den
      if keyword_set(errlum) then masserr1 = errlum else masserr1=0d
      if keyword_set(errden) then masserr2 = errden else masserr2=0d
      masserr = sqrt(masserr1^2d + reverse(masserr2)^2d)
      return, mass
   endelse

end
