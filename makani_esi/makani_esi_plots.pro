pro makani_esi_plots

   bad = 1d99
   
   ; Planck cosmology
   dist = drt_plumdist(0.459d)

   ; label apertures
   naps = 7
   apnames = ['nuc-3"','nuc-2"','1.7"-NE','1.7"-SW','3.3"-W','3.7"-30degNoW',$
      '3.7"-30degSofW']

   ; fluxes, line ratios
   restore,'/Users/drupke/specfits/esi/makani/makani.lin.xdr'
   ; redshifts
   restore,'/Users/drupke/specfits/esi/makani/makani.lininit.xdr'
   ; stellar properties
   restore,'/Users/drupke/specfits/esi/makani/makani.cont.xdr'
   ; Get data
   linelab = 1b
   linelist = ifsf_linelist(/all,/vacuum,linelab=linelab,/quiet)

   ; fluxes and errors of lines
   fx = hash()
   fxe = hash()
   ecfx = hash()
   ecfxe = hash()
   ; gas zs and sigmas
   z = dindgen(naps,2) + bad
   sig = dindgen(naps,2) + bad
   ; populate stellar properties
   stel_z = contcube.stel_z[1:naps]
   stel_sigma = contcube.stel_sigma[1:naps]
   stel_ebv = contcube.stel_ebv[1:naps]
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
   dens_s2 = dindgen(naps,3) + bad
   dens_s2errlo = dindgen(naps,3) + bad
   dens_s2errhi = dindgen(naps,3) + bad
   dens_o2 = dindgen(naps,3) + bad
   dens_o2errlo = dindgen(naps,3) + bad
   dens_o2errhi = dindgen(naps,3) + bad
   ; Halpha luminosities
   loghalum = dindgen(naps,3) + bad
   loghalumerrlo = dindgen(naps,3) + bad
   loghalumerrhi = dindgen(naps,3) + bad
   ; Ionized gas masses, derived from Halpha
   hamass = dindgen(naps,3) + bad
   hamasserrlo = dindgen(naps,3) + bad
   hamasserrhi = dindgen(naps,3) + bad
   ; Halpha star formation rates
   hasfr = dindgen(naps,3) + bad
   hasfrerr = dindgen(naps,3) + bad
   
   ; output text table
   openw,luntxt1,'/Users/drupke/ifs/esi/docs/makani_esi_gasprop.txt',/get_lun
   openw,luntxt2,'/Users/drupke/ifs/esi/docs/makani_esi_gasvel.txt',/get_lun
   openw,luntxt3,'/Users/drupke/ifs/esi/docs/makani_esi_stars.txt',/get_lun
   openw,luntxt4,'/Users/drupke/ifs/esi/docs/makani_esi_chandra.txt',/get_lun
   printf,luntxt1,'Name','ebv','err','logrho_s2','errlo','errhi','logrho_o2',$
      'errlo','errhi','hiimass','errlo','errhi','SFR','err',$
      format='(A-20,13A10)
   printf,luntxt2,'Name','z','siginit',format='(A-20,2A10)
   printf,luntxt3,'Name','z','sigma','ebv',format='(A-20,3A10)
   printf,luntxt4,'Name','[OII]/Ha',format='(A-20,1A10)
   
   for iap=0,naps-1 do begin
      iapp1 = iap+1
      ; populate fluxes and errors of lines
      foreach key,emlflx['fc1'].keys() do begin
         if iap eq 0 then begin
            fx[key] = dindgen(naps,3)+bad
            fxe[key] = dindgen(naps,3)+bad
         endif
         if iap eq 0 or iap eq 1 then begin
            fx[key,iap,*] = [emlflx['fc1',key,iapp1],$
               emlflx['fc2',key,iapp1],$
               emlflx['ftot',key,iapp1]]
            fxe[key,iap,*] = [emlflxerr['fc1',key,iapp1],$
               emlflxerr['fc2',key,iapp1],$
               emlflxerr['ftot',key,iapp1]]
         endif else begin
            fx[key,iap,*] = [emlflx['fc1',key,iapp1],$
               bad,$
               emlflx['fc1',key,iapp1]]
            fxe[key,iap,*] = [emlflxerr['fc1',key,iapp1],$
               bad,$
               emlflxerr['fc1',key,iapp1]]
         endelse
      endforeach
      ; populate gas zs and sigmas
      z_tmp = reform(emlz['[OII]3729',iapp1,*,*])
      z[iap,0:n_elements(z_tmp)-1] = z_tmp
      sig_tmp = reform(emlsiginit['[OII]3729',iapp1,*,*])
      sig[iap,0:n_elements(z_tmp)-1] = sig_tmp
      for i=0,n_elements(z_tmp)-1 do begin
         if z[iap,i] ne bad then begin
            name = apnames[iap]+'_'+string(i+1,format='(I0)')
            printf,luntxt2,name,z[iap,i],sig[iap,i],format='(A-20,D10.6,D10.2)
         endif
      endfor
   endfor

   ; compute extincted line ratios
   lrs = ifsf_lineratios(fx,fxe,linelist,errlo=lrsel,errhi=lrseh,errlin=lrselin)
  
   for iap=0,naps-1 do begin
      iapp1 = iap+1
      ; compute E(B-V)
      foreach key,fx.keys() do begin
         if iap eq 0 then begin
            ecfx[key] = dindgen(naps,3)+bad
            ecfxe[key] = dindgen(naps,3)+bad
         endif
         igd = where(fx[key,iap,*] ne bad AND lrs['ebv',iap,*] ne bad,ctgd)
         if ctgd gt 0 then begin
            ec_tmp = drt_dustcor_ccm(linelist[key],fx[key,iap,igd],$
               lrs['ebv',iap,igd],$
               ebverr=lrselin['ebv',iap,igd],fluxerr=fxe[key,iap,igd])
            ecfx[key,iap,igd] = ec_tmp[*,0]
            ecfxe[key,iap,igd] = ec_tmp[*,1]
         endif
      endforeach
   endfor
   
   ; compute unextincted line ratios
   eclrs = ifsf_lineratios(ecfx,ecfxe,linelist,errlo=eclrsel,errhi=eclrseh,$
      errlin=eclrselin,/lronly)

   drt_voplot,lrs,'/Users/drupke/ifs/esi/plots/makani_esi_vo_noextcor.eps',$
      errlo=lrsel,errhi=lrseh ;,scol=['cg1','cg2','cg3']
   drt_voplot,eclrs,'/Users/drupke/ifs/esi/plots/makani_esi_vo.eps',$
      errlo=eclrsel,errhi=eclrseh ;,scol=['cg1','cg2','cg3']

   igds2 = where(lrs['s2'] ne bad)
   dens_s2[igds2] = drt_nebden(10d^lrs['s2',igds2],'s2',$
      err=lrselin['s2',igds2],denerrlo=dens_s2errlo_tmp,$
      denerrhi=dens_s2errhi_tmp)
   dens_s2errlo[igds2] = dens_s2errlo_tmp
   dens_s2errhi[igds2] = dens_s2errhi_tmp
   igdo2 = where(lrs['o2'] ne bad)
   dens_o2[igdo2] = drt_nebden(10d^lrs['o2',igdo2],'o2',$
      err=lrselin['o2',igdo2],denerrlo=dens_o2errlo_tmp,$
      denerrhi=dens_o2errhi_tmp)
   dens_o2errlo[igdo2] = dens_o2errlo_tmp
   dens_o2errhi[igdo2] = dens_o2errhi_tmp

   ; fluxes, luminosities, etc.
   igdha = where(ecfx['Halpha'] ne bad AND ecfxe['Halpha'] ne bad)
   haflux = ecfx['Halpha',igdha]*1d-16*1d-3
   hafluxerr = ecfxe['Halpha',igdha]*1d-16*1d-3
   halum = drt_linelum(haflux,dist,/ergs,err=hafluxerr,$
      lumerr=halumerr)
   loghalum[igdha] = alog10(halum)
   loghalumerrlo[igdha] = loghalum[igdha] - alog10(halum - halumerr)
   loghalumerrhi[igdha] = alog10(halum+halumerr)-loghalum[igdha]
   hamass[igdha] = drt_hiimass(loghalum[igdha],dens_s2[igdha],$
      errlum=[[loghalumerrlo[igdha]],[loghalumerrhi[igdha]]],$
      errden=[[dens_s2errlo[igdha]],[dens_s2errhi[igdha]]],/log,$
      masserr=hamasserr_tmp)
   hamasserrlo[igdha] = hamasserr_tmp[*,0]
   hamasserrhi[igdha] = hamasserr_tmp[*,1]
      
   ; https://ned.ipac.caltech.edu/level5/Sept12/Calzetti/Calzetti1_2.html
   hasfr[igdha] = 5.5d-42 * halum
   hasfrerr[igdha] = 5.5d-42 * halumerr

   for iap=0,naps-1 do begin
      if iap eq 0 or iap eq 1 then $
         name = apnames[iap]+['_1','_2','_T'] $
      else name = apnames[iap]+['_T']
      for i=0,n_elements(name)-1 do begin
         printf,luntxt1,name[i],lrs['ebv',iap,i],lrsel['ebv',iap,i],$
            dens_s2[iap,i],dens_s2errlo[iap,i],dens_s2errhi[iap,i],$
            dens_o2[iap,i],dens_o2errlo[iap,i],dens_o2errhi[iap,i],$
            hamass[iap,i],hamasserrlo[iap,i],hamasserrhi[iap,i],$
            hasfr[iap,i],hasfrerr[iap,i],$
            format='(A-20,13D10.4)
         if fx['Halpha',iap,i] ne bad AND $
            fx['[OII]3726+[OII]3729',iap,i] ne bad then $
            printf,luntxt4,name[i],$
               fx['[OII]3726+[OII]3729',iap,i]/fx['Halpha',iap,i],$
               format='(A-20,D10.4)'
      endfor
      printf,luntxt3,apnames[iap],stel_z[iap],stel_sigma[iap],stel_ebv[iap],$
         format='(A-20,D10.6,2D10.2)
   endfor

   free_lun,luntxt1
   free_lun,luntxt2
   free_lun,luntxt3
   free_lun,luntxt4

end
