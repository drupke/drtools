; docformat = 'rst'
;
;+
;
; Produces plots showing the entire data spectra, as well as zoomed-in
; portions.
;
; :Categories:
;    IFSFIT
;
; :Returns:
;    Postscript file with plots.
;
; :Params:
;    wave: in, required, type=dblarr(N)
;      Observed wavelengths.
;    flux: in, required, type=dblarr(N)
;      Observed fluxes.
;    zsys: in, required, type=double
;      Galaxy redshift.
;    outfile: in, required, type=str
;      Filename of output plot.
;
; :Keywords:
;    ignorewave: in, optional, type=dblarr(N,2)
;      List of N regions to ignore in determining vertical plot range. Each 
;      region has a lower and upper wavelength limit.
;    title: in, optional, type=string
;      Plot title.
;
; :Author:
;    David S. N. Rupke::
;      Rhodes College
;      Department of Physics
;      2000 N. Parkway
;      Memphis, TN 38104  
;      andto94@gmail.com
;
; :History:
;    ChangeHistory::
;      2017dec18, DSNR, created (copied from cos_specplot.pro)
;
; :Copyright:
;    Copyright (C) 2017 David S. N. Rupke
;
;    This program is free software: you can redistribute it and/or
;    modify it under the terms of the GNU General Public License as
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
PRO drt_specplot, wave, flux, zsys, outfile, ignorewave=ignorewave, title=title

;  List of emission/absorption lines and their corresponding wavelengths.
   linelab = 1b
   lines = ifsf_linelist(!NULL,linelab=linelab,/all)
;  Get lists from hashes
   LineWavelength_list = lines.values()
   LineLabel_list = linelab.values()
;  Convert lists to arrays    
   LineWavelength = LineWavelength_list.toarray()
   LineLabel = LineLabel_list.toarray()
  
;  Shifts the absorption/emission waveelngths by the galaxy's redshift
   ShiftedLines= LineWavelength*(1d + zsys)
  
;  Avoid line label collisions
   closethresh = 2d
   nlines = n_elements(LineWavelength)
   linesall = ShiftedLines
   linelaball = LineLabel
   linecolall = strarr(nlines)+'Blue'
   lineyfracall = dblarr(nlines)+0.05d
   isort_linesall = sort(linesall)
   sort_linesall = linesall[isort_linesall]
   sort_linelaball = linelaball[isort_linesall]
   sort_linecolall = linecolall[isort_linesall]
   sort_lineyfracall = lineyfracall[isort_linesall]
   dlines = sort_linesall[1:nlines-1] - sort_linesall[0:nlines-2]
   iclose = where(dlines lt closethresh,ctclose)
   if ctclose gt 0 then begin
      for i=0,ctclose-1 do begin
         sort_linelaball[iclose[i]] = $
            sort_linelaball[iclose[i]]+', '+sort_linelaball[iclose[i]+1]
         sort_linelaball[iclose[i]+1] = ''
         sort_lineyfracall[iclose[i]] = 0.05d
      endfor
   endif

   relativeflux = flux/median(flux)

;  Get rid of other regions
   if keyword_set(ignorewave) then begin
      signore = size(ignorewave)
      if signore[0] eq 1 then nignore=1 else nignore=signore[2]
      for i=0,nignore-1 do begin
         igd_tmp = where(wave le ignorewave[0,i] OR wave ge ignorewave[1,i])
         if i eq 0 then iignore_not = igd_tmp else $
         iignore_not = cgsetintersection(iignore_not,igd_tmp)
      endfor
      igdwave = iignore_not
   endif else igdwave = indgen(n_elements(wave))

;  Arrays with only good regions
   relativeflux_gd = relativeflux[igdwave]
   wave_gd = wave[igdwave]
  
;  Set the x/y-range of the big spectra
   xran=[min(wave),max(wave)]
   yran=[0,Max(relativeflux_gd)]

;  Acquires the range of the wavelength values
   baserange=Max(wave)-Min(wave)

   xsize=7.5
   ysize=8.5 
   xoffset=(8.5-xsize)/2.0d
   yoffset=(11.0-ysize)/2.0d

;  Various plotting defaults  
   plotchars=0.75
   labchars=0.5
   labchart=1
   x1=0.07
   x2=0.99
  
;  Opens a Postscript value to draw plots on
   cgPS_OPEN,outfile,$
            /Encapsulated,scale_factor=1,Charsize=.5,/NOMATCH,$
            xsize=xsize, ysize=ysize, xoffset=xoffset, yoffset=yoffset, /Inches
   !Y.THICK=2
   !X.THICK=2
   !Y.TICKFORMAT='(F0.1)'
   !Y.MINOR = 2
   !X.TICKLEN=0.05
   !Y.TICKLEN=0.01
    
   if keyword_set(title) then $
      cgText,.5,.98,Title,alignment=.5, Charsize = 1,/norm
  
;  Full range plot
   cgplot, wave, relativeflux, xstyle=1, ystyle=1, yran=yran,xran=xran,$
           xtit='Observed Wavelength ($\Angstrom$)',$
           ytit='F!I$\lambda$!N/10!E-14!N (ergs s!E-1 !Ncm!E-2 !N$\Angstrom$!E-1!N)', $
           Position = [x1,.86,x2,.97], CHARSIZE=plotchars,thick=0.5,/NoErase,$
           xticklen=0.05,yticklen=0.01
   index=Where(sort_linesall lt xran[1] AND sort_linesall gt xran[0],ctgd)
   if ctgd gt 0 then begin
      gdline = sort_linesall[index]
      gdlinelab = sort_linelaball[index]
      gdlinecol = sort_linecolall[index]
      gdlineyfrac = sort_lineyfracall[index]
      FOR M = 0,ctgd-1 DO BEGIN
         cgoplot,[gdline[M],gdline[M]],yran, color = gdlinecol[m], thick=1
;         cgTEXT,gdline[m]-.1,yran[1]-(yran[0]+yran[1])*gdlineyfrac[m],$
;                gdlinelab[m],/Data,ORIENT=90d,CHARSIZE = labchars,$
;                charthick=labchart,align=1
      endfor
   endif


;  Plots 3 zoomed-in regions, one after another, $
;  along with absorption/emission lines and labels

   y1 = [0.56,0.29,0.02]
   y2 = [0.805,0.535,0.265]
   for i=0,2 do begin
      xran = [Min(wave)+(double(i)/3)*baserange,$
              Min(wave)+(double(i+1)/3)*baserange]
      fluxtmp = relativeflux_gd[value_locate(wave_gd,xran[0]):$
                                value_locate(wave_gd,xran[1])]
      min = min(fluxtmp)
      max = max(fluxtmp)
      fluxran = max - min
      yran = [min - 0.05d*fluxran, max+0.05d*fluxran]
      cgplot, wave, relativeflux, xstyle=1, ystyle=1, $
              xran=xran, axiscolor='Black',color='Black',yran=yran, $
              Position = [x1,y1[i],x2,y2[i]], /NoErase, CHARSIZE=plotchars,$
              thick=0.5
      index=Where(sort_linesall lt xran[1] AND sort_linesall gt xran[0],ctgd)
      if ctgd gt 0 then begin
         gdline = sort_linesall[index]
         gdlinelab = sort_linelaball[index]
         gdlinecol = sort_linecolall[index]
         gdlineyfrac = sort_lineyfracall[index]
         FOR M = 0,ctgd-1 DO BEGIN
            cgoplot,[gdline[M],gdline[M]],yran, color = gdlinecol[m], thick=1
            cgTEXT,gdline[m]-.1,yran[1]-(yran[0]+yran[1])*gdlineyfrac[m],$
                   gdlinelab[m],/Data,ORIENT=90d,CHARSIZE = labchars,$
                   charthick=labchart,align=1
         endfor
      endif
         
   endfor

   CGPS_CLOSE

END
