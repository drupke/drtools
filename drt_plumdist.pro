; docformat = 'rst'
;
;+
;
; Compute luminosity and angular size distances using Planck 2018 cosmology.
;
; :Categories:
;    none
;
; :Returns:
;    Luminosity distance in Mpc.
;
; :Params:
;    z: in, required, type=double
;
; :Keywords:
;    angdist: out, optional, type=double
;      Also compute angular size distance, put in this variable
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
;      2022jan28, DSNR, created
;    
; :Copyright:
;    Copyright (C) 2022 David S. N. Rupke
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
function drt_plumdist,z,angdist=angdist

   dl = lumdist(z,H0=67.4d,Omega_M=0.315d,Lambda0=0.685d)
   if keyword_set(angdist) then angdist = dl/(1d + z)^2d

   return,dl
end
