; docformat = 'rst'
;
;+
;
; Compute intrinsic flux and its error from the Cardelli, Clayton, &
; Mathis (1989) extinction curve.
;
; Root equation is
;     m_i = m_o - (A/A_V * E(B-V) * R_V)
;
; Substitute fluxes for mags. to get the equations below.
;
; :Categories:
;    DRTOOLS
;
; :Returns:
;    If neither FLUXERR nor EBVERR is input, then corrected fluxes as
;    dblarr(N), where N is the number of observed fluxes input. If
;    FLUXERR or EBVERR is input, then dblarr(N,2), where the first row
;    is corrected flux and the second is its error.
;
; :Params:
;    lambda: in, required, type=double or dblarr(2)
;        If one number is input, wavelength at which to compute the
;        intrinsic flux. If two, and the keyword RELATIVE is set, then
;        the relative correction between two wavelengths is applied.
;    flux: in, required, type=dblarr(N)
;        Observed fluxes to correct.
;    ebv: in, required, type=dblarr(N)
;        Observed values of E(B-V), one per flux.
;
; :Keywords:
;    log: in, optional, type=boolean
;        If set, input/output fluxes and errors are one-sided in the log.
;    rv: in, optional, type=double, default=3.1
;        Total selective extinction A(V)/E(B-V).
;    fluxerr: in, optional, type=dblarr(N)
;        Errors in input fluxes. If set, corrected flux errors are
;        calculated and output.
;    ebverr: in, optional, type=dblarr(2) or dblarr(N,2)
;        Errors in E(B-V). If set, corrected flux errors are
;        calculated and output.
;    tran: in, optional, type=boolean
;        If this keyword is set, there is only one flux, and errors
;        are computed, transpose the output from dblarr(1,2) to
;        dblarr(2).
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
;      2009sep08, DSNR, created
;      2021dec17, DSNR, update E(B-V) error calculation
;      2023oct12, DSNR, check that E(B-V) != bad
;
; :Copyright:
;    Copyright (C) 2009-2023 David S. N. Rupke
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
function drt_dustcor_ccm,lambda,flux,ebv,rv=rv,fluxerr=fluxerr,$
                         ebverr=ebverr,relative=relative,tran=tran,$
                         log=log

  bad = 1d99

  if ~ keyword_set(rv) then rv=3.1d

  alamav = drt_extcurve_ccm(lambda,rv=rv)
  if keyword_set(relative) then coeff = 0.4d *(alamav[0]-alamav[1])*rv $
  else coeff = 0.4d * alamav[0] * rv

; Only apply extinction correction if calculated extinction is greater
; than 0.
  posebv = where(ebv ge 0 AND ebv ne bad,ctpos)
  negebv = where(ebv lt 0,ctneg)

; Fluxes
  fint = dblarr(n_elements(ebv))
  finterr = dblarr(n_elements(ebv))
  if ctpos gt 0 then $
     if not keyword_set(log) then $
        fint[posebv] = flux[posebv] * 10d^(coeff * ebv[posebv]) $
     else $
        fint[posebv] = flux[posebv] + coeff*ebv[posebv]
  if ctneg gt 0 then $
     fint[negebv] = flux[negebv]

; Errors
  doerr=0
  if keyword_set(fluxerr) then begin
     doerr=1
     if ctpos gt 0 then $
        if not keyword_set(log) then $
           finterr[posebv] += fluxerr[posebv] * 10d^(coeff * ebv[posebv]) $
        else $
           finterr[posebv] += fluxerr[posebv]
     if ctneg gt 0 then $
        finterr[negebv] = fluxerr[negebv]
  endif
  if keyword_set(ebverr) then begin
     doerr=1
     if ctpos gt 0 then begin
;       New approach and old approach
;       First one is from dfint = sqrt[(dfint/df*finterr)^2+(dfint/debv*ebverr)^2)]
;       second one is ... ?
        if not keyword_set(log) then begin
           finterr[posebv] = finterr[posebv]^2d
           finterr[posebv] += (coeff * fint[posebv] * ebverr[posebv])^2d
           finterr[posebv] = sqrt(finterr[posebv])
;        finterr[posebv] += ebverr[posebv] * flux[posebv] * coeff * $
;                           alog(10d) * exp(coeff * alog(10d) * $
;                                           ebv[posebv])
        endif else begin
           finterr[posebv] = sqrt(finterr[posebv]^2d + (coeff*ebverr[posebv])^2d)
        endelse
     endif
     if ctneg gt 0 then $
        finterr[negebv] = fluxerr[negebv]
  endif
  if doerr then begin
     fint = [[fint],[finterr]]
     if n_elements(flux) eq 1 AND keyword_set(tran) then fint = transpose(fint)
  endif

  return,fint

end
