; docformat = 'rst'
;
;+
;
; Plot data on Veilleux/Osterbrock 1987 (aka BPT) plots.
;
; :Categories:
;    DRTOOLS
;
; :Returns:
;    Postscript plots.
;
; :Params:
;    lrat: in, required, type=hash
;    outfile: in, required, type=string
;
; :Keywords:
;    errlo: in, optional, type=hash
;    errhi: in, optional, type=hash
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
;      2009xxxYY  DSNR  created
;      2010may18  DSNR  added error representation
;      2012may08  DSNR  rewrote using CG
;      2019apr10  DSNR  a few bugfixes
;      2021dec16  DSNR  ported to DRT library
;
; :Copyright:
;    Copyright (C) 2009--2021 David S. N. Rupke
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
pro drt_voplot,lrat,fileroot,errlo=errlo,errhi=errhi,scol=scol

   cgps_open,fileroot,xsize=7.5,ysize=7.5,/inches,/encap,$
      charsize=1,default_thick=1,/qui,/nomatch

   if keyword_set(errlo) and keyword_set(errhi) then doerr=1b else doerr=0b

   types = lrat.keys()
   ssize = dblarr(n_elements(lrat[types[0]]))+1d
   stype = intarr(n_elements(lrat[types[0]]))+16
   if ~ keyword_set(scol) then scol = strarr(n_elements(types[0]))+'Black'

; [NII]/Halpha vs. [OIII]/Hb

   xlab = 'n2ha'
   ylab = 'o3hb'
   if lrat.haskey(xlab) AND lrat.haskey(ylab) then begin

      xtit='log([NII]/H$\alpha$)'
      ytit='log([OIII]/H$\beta$)'
      xran=[-1.99d,0.99d]
      yran=[-1.19d,1.49d]
      cgplot,[0],xran=xran,yran=yran,ytit=ytit,xtit=xtit,$
         /xsty,/ysty,layout=[2,2,1],/nodata
      ;  Kauffmann et al. 2003 dividing lines
      kf_n2ha1 = 0.05d*dindgen(100)-5d
      kf_o3hb1 = 1.3d + 0.61d / (kf_n2ha1 - 0.05d)
      cgoplot,kf_n2ha1,kf_o3hb1,linestyle=2
      ; Kewley et al. 2001/6 dividing lines
      xkew1 = 0.05d*indgen(110)-5d
      ykew1 = 0.61d / (xkew1-0.47d)+1.19d
      cgoplot,xkew1,ykew1,linestyle=1
      ; data
      x = lrat[xlab]
      y = lrat[ylab]
      for i=0,n_elements(x)-1 do begin
         if doerr then begin
            xel = errlo[xlab]
            xeh = errhi[xlab]
            yel = errlo[ylab]
            yeh = errhi[ylab]
            cgoplot,x[i],y[i],psym=stype[i],symsize=ssize[i],color=scol[i],$
               err_xlo=xel[i],err_xhi=xeh[i],err_ylo=yel[i],err_yhi=yeh[i],$
               err_color='Grey',err_thick=2,/err_clip,err_width=0d
         endif else begin
            cgoplot,x[i],y[i],psym=stype[i],symsize=ssize[i],color=scol[i]
         endelse
      endfor

   endif

; [SII]/Halpha vs. [OIII]/Hb

   xlab = 's2ha'
   ylab = 'o3hb'
   if lrat.haskey(xlab) AND lrat.haskey(ylab) then begin

      xtit='log([SII]/H$\alpha$)'
      ytit='log([OIII]/H$\beta$)'
      xran=[-1.19d,0.79d]
      yran=[-1.19d,1.49d]
      cgplot,[0],xran=xran,yran=yran,ytit=ytit,xtit=xtit,$
         /xsty,/ysty,layout=[2,2,2],/nodata
      xkew1 = 0.05*indgen(105)-5
      ykew1 = 0.72d / (xkew1-0.32d)+1.30d
      xkew2 = 0.5d*indgen(2)-0.4d
      ykew2 = 1.89d*xkew2+0.76d
      cgoplot,xkew1,ykew1,linestyle=1
      cgoplot,xkew2,ykew2,linestyle=1
      x = lrat[xlab]
      y = lrat[ylab]
      for i=0,n_elements(x)-1 do begin
         if doerr then begin
            xel = errlo[xlab]
            xeh = errhi[xlab]
            yel = errlo[ylab]
            yeh = errhi[ylab]
            cgoplot,x[i],y[i],psym=stype[i],symsize=ssize[i],color=scol[i],$
               err_xlo=xel[i],err_xhi=xeh[i],err_ylo=yel[i],err_yhi=yeh[i],$
               err_color='Grey',err_thick=2,/err_clip,err_width=0d
         endif else begin
            cgoplot,x[i],y[i],psym=stype[i],symsize=ssize[i],color=scol[i]
         endelse
      endfor
  
   endif

; [OI]/Halpha vs. [OIII]/Hb

   xlab = 'o1ha'
   ylab = 'o3hb'
   if lrat.haskey(xlab) AND lrat.haskey(ylab) then begin

      xtit='log([OI]/H$\alpha$)'
      ytit='log([OIII]/H$\beta$)'
      xran=[-2.19d,-0.01d]
      yran=[-1.19d, 1.49d]
      cgplot,[0],xran=xran,yran=yran,ytit=ytit,xtit=xtit,/xsty,/ysty,$
         layout=[2,2,3],/nodata
      xkew1 = 0.05*indgen(85)-5
      ykew1 = 0.73d / (xkew1+0.59d)+1.33d
      xkew2 = 0.5d*indgen(2)-1.1d
      ykew2 = 1.18d*xkew2 + 1.30d
      cgoplot,xkew1,ykew1,linestyle=1
      cgoplot,xkew2,ykew2,linestyle=1
      x = lrat[xlab]
      y = lrat[ylab]
      for i=0,n_elements(x)-1 do begin
         if doerr then begin
            xel = errlo[xlab]
            xeh = errhi[xlab]
            yel = errlo[ylab]
            yeh = errhi[ylab]
            cgoplot,x[i],y[i],psym=stype[i],symsize=ssize[i],color=scol[i],$
               err_xlo=xel[i],err_xhi=xeh[i],err_ylo=yel[i],err_yhi=yeh[i],$
               err_color='Grey',err_thick=2,/err_clip,err_width=0d
         endif else begin
            cgoplot,x[i],y[i],psym=stype[i],symsize=ssize[i],color=scol[i]
         endelse
      endfor

   endif

; [OI]/Halpha vs. [OIII]/[OII]

   xlab = 'o1ha'
   ylab = 'o3o2'
   if lrat.haskey(xlab) AND lrat.haskey(ylab) then begin

      xtit='log([OI]/H$\alpha$)'
      ytit='log([OIII]/[OII])'
      xran=[-2.39d,-0.01d]
      yran=[-1.49d,0.99d]
      cgplot,[0],xran=xran,yran=yran,ytit=ytit,xtit=xtit,$
         /xsty,/ysty,layout=[2,2,4],/nodata
      kw_o1ha_a = 0.05*indgen(100)-2.5
      kw_o3o2_a = -1.701*kw_o1ha_a-2.163
      kw_o1ha_b = 0.05*indgen(100)-1.05
      kw_o3o2_b = kw_o1ha_b+0.7
      cgoplot,kw_o1ha_a,kw_o3o2_a,linestyle=1
      cgoplot,kw_o1ha_b,kw_o3o2_b,linestyle=1
      x = lrat[xlab]
      y = lrat[ylab]
      for i=0,n_elements(x)-1 do begin
         if doerr then begin
            xel = errlo[xlab]
            xeh = errhi[xlab]
            yel = errlo[ylab]
            yeh = errhi[ylab]
            cgoplot,x[i],y[i],psym=stype[i],symsize=ssize[i],color=scol[i],$
               err_xlo=xel[i],err_xhi=xeh[i],err_ylo=yel[i],err_yhi=yeh[i],$
               err_color='Grey',err_thick=2,/err_clip,err_width=0d
         endif else begin
            cgoplot,x[i],y[i],psym=stype[i],symsize=ssize[i],color=scol[i]
         endelse
      endfor

   endif

   cgps_close

end
