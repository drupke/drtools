pro orcs_plots

   bad = 1d99
   
   speryr = 24d*3600d*365.25d    ; seconds in a year
   lsun = 3.826d33               ; solar luminosities, erg/s
   msun = 1.989e33               ; solar mass, g
   c_cms = 2.99792d10

   plotdir = '/Users/drupke/ifs/gmos/plots/orcs/'
   plotquantum = 2.5 ; in inches


   ; Planck cosmology
   angdist=1b
   dist = drt_plumdist(0.451d,angdist=angdist)

   kpc_per_as = angdist*1000d/206265d
   gmosps = 0.1614d
   kpc_per_pix = gmosps * kpc_per_as

   ; fluxes, line ratios
   restore,'/Users/drupke/specfits/gmos/orc4c-ql/orc4c-ql.lin.xdr'
   ; redshifts
   restore,'/Users/drupke/specfits/gmos/orc4c-ql/orc4c-ql.lininit.xdr'
   ; stellar properties
   restore,'/Users/drupke/specfits/gmos/orc4c-ql/orc4c-ql.cont.xdr'
   ; apertures
   naps = 1
   ; Get data
   linelab = 1b
   linelist = ifsf_linelist(/all,linelab=linelab,/quiet)

   ; Compute CVDF velocities
   emlcvdfvel = ifsf_cmpcvdfvals(emlcvdf,emlflx,emlflxerr)

   ; fluxes and errors of lines
   fx = hash()
   fxe = hash()
   ecfx = hash()
   ecfxe = hash()
   ; gas zs and sigmas
   ;z = dindgen(naps,2) + bad
   ;sig = dindgen(naps,2) + bad
   ; populate stellar properties
   ;stel_z = contcube.stel_z[apoff:apoff+naps-1]
   ;stel_sigma = contcube.stel_sigma[apoff:apoff+naps-1]
   ;stel_ebv = contcube.stel_ebv[apoff:apoff+naps-1]
   ; line ratios
   lrs = hash()
   lrsel = hash()
   lrseh = hash()
   lrselin = hash()
   eclrs = hash()
   eclrsel = hash()
   eclrseh = hash()
   eclrselin = hash()
   ; densities
   dens_s2 = dindgen(naps) + bad
   dens_s2errlo = dindgen(naps) + bad
   dens_s2errhi = dindgen(naps) + bad
   dens_o2 = dindgen(naps) + bad
   dens_o2errlo = dindgen(naps) + bad
   dens_o2errhi = dindgen(naps) + bad
   ; Halpha luminosities
   loghalum = dindgen(naps) + bad
   loghalumerrlo = dindgen(naps) + bad
   loghalumerrhi = dindgen(naps) + bad
   ; Ionized gas masses, derived from Halpha
   ;hamass = dindgen(naps) + bad
   ;hamasserrlo = dindgen(naps) + bad
   ;hamasserrhi = dindgen(naps) + bad
   ; Halpha star formation rates
   ;hasfr = dindgen(naps) + bad
   ;hasfrerr = dindgen(naps) + bad
   
;   ; output text table
;   openw,luntxt1,'/Users/drupke/ifs/esi/docs/makani_esi_gasprop.txt',/get_lun
;   openw,luntxt2,'/Users/drupke/ifs/esi/docs/makani_esi_gasvel.txt',/get_lun
;   openw,luntxt3,'/Users/drupke/ifs/esi/docs/makani_esi_stars.txt',/get_lun
;   openw,luntxt4,'/Users/drupke/ifs/esi/docs/makani_esi_chandra.txt',/get_lun
;   printf,luntxt1,'Ap#','Name','ebv','err','logrho_s2','errlo','errhi','logrho_o2',$
;      'errlo','errhi','hiimass','errlo','errhi','SFR','err',$
;      format='(A-3,A-15,13A10)
;   printf,luntxt2,'Name','z','siginit',format='(A-20,2A10)
;   printf,luntxt3,'Name','z','sigma','ebv',format='(A-20,3A10)
;   printf,luntxt4,'Name','[OII]/Ha',format='(A-20,1A10)
   
   for iap=0,naps-1 do begin
      ncomp=1
      ; populate fluxes and errors of lines
      foreach key,emlflx['fc1'].keys() do begin
         if iap eq 0 then begin
            fx[key] = dindgen(naps)+bad
            fxe[key] = dindgen(naps)+bad
         endif
         fx[key,iap] = emlflx['fc1',key]
         fxe[key,iap] = emlflxerr['fc1',key]
      endforeach
      ; populate gas zs and sigmas
;      z_tmp = reform(emlz['[OII]3729',iapp1,*,*])
;      z[iap,0:n_elements(z_tmp)-1] = z_tmp
;      sig_tmp = reform(emlsiginit['[OII]3729',iapp1,*,*])
;      sig[iap,0:n_elements(z_tmp)-1] = sig_tmp
;      for i=0,n_elements(z_tmp)-1 do begin
;         if z[iap,i] ne bad then begin
;            name = apnames[iap]+'_'+string(i+1,format='(I0)')
;            printf,luntxt2,name,z[iap,i],sig[iap,i],format='(A-20,D10.6,D10.2)
;         endif
;      endfor
   endfor

   ; compute extincted line ratios
   lrs = ifsf_lineratios(fx,fxe,linelist,errlo=lrsel,errhi=lrseh,errlin=lrselin)
  
   for iap=0,naps-1 do begin
      ; compute E(B-V)
      foreach key,fx.keys() do begin
         if iap eq 0 then begin
            ecfx[key] = dindgen(naps)+bad
            ecfxe[key] = dindgen(naps)+bad
         endif
         igd = where(fx[key,iap] ne bad AND lrs['ebv',iap] ne bad,ctgd)
         if ctgd gt 0 then begin
            ec_tmp = drt_dustcor_ccm(linelist[key],fx[key,iap],$
               lrs['ebv',iap],ebverr=lrselin['ebv',iap],$
               fluxerr=fxe[key,iap])
            ecfx[key,iap] = ec_tmp[*,0]
            ecfxe[key,iap] = ec_tmp[*,1]
         endif
      endforeach
   endfor
   
   print,'E(B-V)'
   print,string(lrs['ebv'],format='(D0.2)')+'+/-'+string(lrselin['ebv'],format='(D0.2)')
   
   ; compute unextincted line ratios
   eclrs = ifsf_lineratios(ecfx,ecfxe,linelist,errlo=eclrsel,errhi=eclrseh,$
      errlin=eclrselin,/lronly)

   drt_voplot,lrs,plotdir+'orcs_vo_noextcor.eps',errlo=lrsel,errhi=lrseh
   drt_voplot,eclrs,plotdir+'orcs_vo.eps',errlo=eclrsel,errhi=eclrseh
   
   ; e- density
   dens_o2 = dindgen(naps) + bad
   dens_o2errlo = dindgen(naps) + bad
   dens_o2errhi = dindgen(naps) + bad
   igdo2 = where(lrs['o2'] ne bad)
   dens_o2[igdo2] = drt_nebden(10d^lrs['o2',igdo2],'o2',$
      err=lrselin['o2',igdo2],denerrlo=dens_o2errlo_tmp,$
      denerrhi=dens_o2errhi_tmp)
   dens_o2errlo[igdo2] = dens_o2errlo_tmp
   dens_o2errhi[igdo2] = dens_o2errhi_tmp
   print,'log [OII] den/cm^-2'
   print,string(dens_o2,format='(D0.2)')+'+'+$
      string(dens_o2errlo,format='(D0.2)')+'-'+$
      string(dens_o2errhi,format='(D0.2)')
   
   ; Halpha luminosities
   loghalum = dindgen(naps) + bad
   loghalumerrlo = dindgen(naps) + bad
   loghalumerrhi = dindgen(naps) + bad
   igdha = where(ecfx['Halpha'] ne bad AND ecfxe['Halpha'] ne bad)
   haflux = ecfx['Halpha',igdha]*1d-17*1d-3 ; fluxes in W m^-2
   hafluxerr = ecfxe['Halpha',igdha]*1d-17*1d-3 ; fluxes in W m^-2
   halum = drt_linelum(haflux,dist,/ergs,err=hafluxerr,$
      lumerr=halumerr)
   loghalum[igdha] = alog10(halum)
   loghalumerrlo[igdha] = loghalum[igdha] - alog10(halum - halumerr)
   loghalumerrhi[igdha] = alog10(halum+halumerr)-loghalum[igdha]
   print,'log L(Ha)/erg/s'
   print,string(loghalum,format='(D0.2)')+'+'+$
      string(loghalumerrlo,format='(D0.2)')+'-'+$
      string(loghalumerrhi,format='(D0.2)')

   ; Ionized gas masses, derived from Halpha
   hamass = dindgen(naps) + bad
   hamasserrlo = dindgen(naps) + bad
   hamasserrhi = dindgen(naps) + bad
   hamass[igdha] = drt_hiimass(loghalum[igdha],2.3d,$
      errlum=[[loghalumerrlo[igdha]],[loghalumerrhi[igdha]]],/log,$
      masserr=hamasserr_tmp)
   hamasserrlo[igdha] = hamasserr_tmp[*,0]
   hamasserrhi[igdha] = hamasserr_tmp[*,1]
   print,'log M(Ha)/Msun (assuming n_e = 200 cm^-2)'
   print,string(hamass,format='(D0.2)')+'+'+$
      string(hamasserrlo,format='(D0.2)')+'-'+$
      string(hamasserrhi,format='(D0.2)')

   ; o3 luminosities
   logo3lum = dindgen(naps) + bad
   logo3lumerrlo = dindgen(naps) + bad
   logo3lumerrhi = dindgen(naps) + bad
   igdo3 = where(ecfx['[OIII]5007'] ne bad AND ecfxe['[OIII]5007'] ne bad)
   o3flux = ecfx['[OIII]5007',igdo3]*1d-17*1d-3 ; fluxes in W m^-2
   o3fluxerr = ecfxe['[OIII]5007',igdo3]*1d-17*1d-3 ; fluxes in W m^-2
   o3lum = drt_linelum(o3flux,dist,/ergs,err=o3fluxerr,$
      lumerr=o3lumerr)
   logo3lum[igdo3] = alog10(o3lum)
   logo3lumerrlo[igdo3] = logo3lum[igdo3] - alog10(o3lum - o3lumerr)
   logo3lumerrhi[igdo3] = alog10(o3lum+o3lumerr)-logo3lum[igdo3]
   print,'log L([OIII])/erg/s'
   print,string(logo3lum,format='(D0.2)')+'+'+$
      string(logo3lumerrlo,format='(D0.2)')+'-'+$
      string(logo3lumerrhi,format='(D0.2)')

   ; Compute intrinsic Hbeta flux across pixel. Divide by arcsec^2 and by cm^-2/as^-2
   ;loghblumsb[igdha] = alog10(halum) - alog10(aplenpixarr[igdha]*esips) - $
   ;   2d*alog10(kpc_per_as * cm_per_pc * 1d3) - alog10(2.86d)
   ;print,loghblumsb
   
;   ; --------------
;   ; Flux table in LaTeX
;   ; --------------
;
;   openw,tab3,'/Users/drupke/ifs/gmos/docs/orcs/orc4c_fluxtab.tex',/get_lun
;   amp = ' & '
;   dslash = ' \\'
;   tablines = ['[OII]3726','[OII]3729','[OII]3726+[OII]3729',$
;      'Hbeta','[OIII]5007','[OI]6300','[NII]6583']
;   ; if set to same as sigcut in IFSFA, just duplicates that and has no effect
;   sigcut = 2d
;   ; errors in log here are from d(log f1 - log f2) = sqrt((df1/f1)^2d+(df2/f2)^2)/ln 10
;   ; d(log x) = 1/(x log 10), and d(g(f1,f2)) = sqrt((dg/df1 * df1)^2d + ...)
;   ; observed fluxes
;   for i=0,n_elements(tablines)-1 do begin
;      for j=1,8 do begin
;         if fx[tablines[i],j,2] ne bad AND $
;            fx[tablines[i],j,2] ne 0d AND $
;            fx['Halpha',j,2] ne bad AND $
;            fx['Halpha',j,2] ne 0 then begin
;            if fx[tablines[i],j,2]/fxe[tablines[i],j,2] gt sigcut then begin
;               drat = sqrt((fxe[tablines[i],j,2]/fx[tablines[i],j,2])^2d + $
;                  (fxe['Halpha',j,2]/fx['Halpha',j,2])^2d) / alog(10d)
;               tabstr += string($
;                  amp,alog10(fx[tablines[i],j,2]/fx['Halpha',j,2]),$
;                  '$\pm$',drat,format='(A3,D7.2,A5,D4.2)')
;            endif else begin
;               tabstr += string(amp,'\nodata',format='(A3,A16)')
;            endelse
;         endif else begin
;            tabstr += string(amp,'\nodata',format='(A3,A16)')
;         endelse
;      endfor
;      if i ne n_elements(tablines)-1 then tabstr += dslash+'\relax' $
;      else tabstr += dslash
;      printf,tab3,tabstr
;   endfor
;   printf,tab3,'\hline'
;   drat1 = fxe['Halpha',0,0]/fx['Halpha',0,0] / alog(10d)
;   drat2 = fxe['Halpha',0,1]/fx['Halpha',0,1] / alog(10d)
;   drat = fxe['Halpha',0,2]/fx['Halpha',0,2] / alog(10d)
;   tabstr = string('log[f(Halpha)]',amp,$
;      alog10(fx['Halpha',0,0])-16d,'$\pm$',drat1,amp,$
;      alog10(fx['Halpha',0,1])-16d,'$\pm$',drat2,amp,$
;      alog10(fx['Halpha',0,2])-16d,'$\pm$',drat,$
;      format='(A30,A3,D7.3,A5,D5.3,A3,D7.3,A5,D5.3,A3,D7.3,A5,D5.3)')
;   for j=1,8 do begin
;      drat = fxe['Halpha',j,2]/fx['Halpha',j,2] / alog(10d)
;      tabstr += string(amp,alog10(fx['Halpha',j,2])-16d,'$\pm$',drat,$
;         format='(A3,D6.2,A5,D4.2)')
;   endfor
;   tabstr += dslash
;   printf,tab3,tabstr
;   printf,tab3,'\hline'
;
;   ; intrinsic fluxes
;   for i=0,n_elements(tablines)-1 do begin
;      drat1 = sqrt((ecfxe[tablines[i],0,0]/ecfx[tablines[i],0,0])^2d + $
;         (ecfxe['Halpha',0,0]/ecfx['Halpha',0,0])^2d) / alog(10d)
;      drat2 = sqrt((ecfxe[tablines[i],0,1]/ecfx[tablines[i],0,1])^2d + $
;         (ecfxe['Halpha',0,1]/ecfx['Halpha',0,1])^2d) / alog(10d)
;      drat = sqrt((ecfxe[tablines[i],0,2]/ecfx[tablines[i],0,2])^2d + $
;         (ecfxe['Halpha',0,2]/ecfx['Halpha',0,2])^2d) / alog(10d)
;      if fx[tablines[i],0,0]/fxe[tablines[i],0,0] gt sigcut then $
;         c1str = string(alog10(ecfx[tablines[i],0,0]/ecfx['Halpha',0,0]),$
;         '$\pm$',drat1,format='(D7.2,A5,D4.2)') $
;      else $
;         c1str = '\nodata'
;      if fx[tablines[i],0,1]/fxe[tablines[i],0,1] gt sigcut then $
;         c2str = string(alog10(ecfx[tablines[i],0,1]/ecfx['Halpha',0,1]),$
;         '$\pm$',drat2,format='(D7.2,A5,D4.2)') $
;      else $
;         c2str = '\nodata'
;      if fx[tablines[i],0,2]/fxe[tablines[i],0,2] gt sigcut then $
;         ctotstr = string(alog10(ecfx[tablines[i],0,2]/ecfx['Halpha',0,2]),$
;         '$\pm$',drat,format='(D7.2,A5,D4.2)') $
;      else $
;         ctotstr = '\nodata'
;      tabstr = string(tablines[i],amp,c1str,amp,c2str,amp,ctotstr,$
;         format='(A20,A3,A16,A3,A16,A3,A16)')
;      for j=1,8 do begin
;         if ecfx[tablines[i],j,2] ne bad AND $
;            ecfx[tablines[i],j,2] ne 0d AND $
;            ecfx['Halpha',j,2] ne bad AND $
;            ecfx['Halpha',j,2] ne 0 then begin
;            if fx[tablines[i],j,2]/fxe[tablines[i],j,2] gt sigcut then begin
;               drat = sqrt((ecfxe[tablines[i],j,2]/ecfx[tablines[i],j,2])^2d + $
;                  (ecfxe['Halpha',j,2]/ecfx['Halpha',j,2])^2d) / alog(10d)
;               tabstr += string($
;                  amp,alog10(ecfx[tablines[i],j,2]/ecfx['Halpha',j,2]),$
;                  '$\pm$',drat,format='(A3,D7.2,A5,D4.2)')
;            endif else begin
;               tabstr += string(amp,'\nodata',format='(A3,A16)')
;            endelse
;         endif else begin
;            tabstr += string(amp,'\nodata',format='(A3,A16)')
;         endelse
;      endfor
;      if i ne n_elements(tablines)-1 then tabstr += dslash+'\relax' $
;      else tabstr += dslash
;      printf,tab3,tabstr
;   endfor
;   printf,tab3,'\hline'
;   tabstr = string('log[L(Halpha)]',amp,$
;      loghalum[0,0],'$^{+',loghalumerrhi[0,0],'}_{-',loghalumerrlo[0,0],'}$',amp,$
;      loghalum[0,1],'$^{+',loghalumerrhi[0,1],'}_{-',loghalumerrlo[0,1],'}$',amp,$
;      loghalum[0,2],'$^{+',loghalumerrhi[0,2],'}_{-',loghalumerrlo[0,2],'}$',$
;      format='(A20,A3,D5.2,A4,D4.2,A4,D4.2,A2,A3,D5.2,A4,D4.2,A4,D4.2,A2,A3,D5.2,A4,D4.2,A4,D4.2,A2)')
;   for j=1,8 do begin
;      tabstr += $
;         string(amp,loghalum[j,2],'$^{+',loghalumerrhi[j,2],'}_{-',loghalumerrlo[j,2],'}$',$
;         format='(A3,D5.2,A4,D4.2,A4,D4.2,A2)')
;   endfor
;   ;tabstr += dslash+'\relax'
;   printf,tab3,tabstr
;
;
;   free_lun,tab3

   print,ecfx
   
end