pro makani_esi_plots

   iap = 1
   apname = '3"nuc'
   iap = 2
   apname = '2"nuc'
   iap = 3
   apname = '+1.7"'
   iap = 4
   apname = '-1.7"'
   iap = 5
   apname = '3"W_cent'
   iap = 7
   apname = '3"W_cent-1.7"'

   aplabp1 = string(iap+1,format='(I0)')
   aplab = string(iap,format='(I0)')

   linelab = 1b
   linelist = ifsf_linelist(/all,/vacuum,linelab=linelab,/quiet)
   restore,'/Users/drupke/specfits/esi/makani/makani_000'+aplabp1+'.xdr'
   restore,'/Users/drupke/specfits/esi/makani/makani.lin.xdr'
   restore,'/Users/drupke/specfits/esi/makani/makani.lininit.xdr'
   restore,'/Users/drupke/specfits/esi/makani/makani.cont.xdr'
   fx = hash()
   fxe = hash()
   foreach key,emlflx['fc1'].keys() do begin
      if iap eq 1 or iap eq 2 then begin
         fx[key] = [emlflx['fc1',key,iap],emlflx['fc2',key,iap],emlflx['ftot',key,iap]]
         fxe[key] = [emlflxerr['fc1',key,iap],emlflxerr['fc2',key,iap],emlflxerr['ftot',key,iap]]
      endif else begin
         fx[key] = [emlflx['fc1',key,iap]]
         fxe[key] = [emlflxerr['fc1',key,iap]]
      endelse
   endforeach
   z = reform(emlz['[OII]3729',iap,*,*])
   sig = reform(emlsiginit['[OII]3729',iap,*,*])
   stel_z = contcube.stel_z[iap]
   stel_sigma = contcube.stel_sigma[iap]
   stel_ebv = contcube.stel_ebv[iap]

   lrs = ifsf_lineratios(fx,fxe,linelist,errlo=lrsel,errhi=lrseh,errlin=lrselin)
   ecfx = hash()
   ecfxe = hash()
   foreach key,fx.keys() do begin
      ec_tmp = drt_dustcor_ccm(linelist[key],fx[key],lrs['ebv'],$
         ebverr=lrselin['ebv'],fluxerr=fxe[key])
      ecfx[key] = ec_tmp[*,0]
      ecfxe[key] = ec_tmp[*,1]
   endforeach
   eclrs = ifsf_lineratios(ecfx,ecfxe,linelist,errlo=eclrsel,errhi=eclrseh,$
      errlin=eclrselin,/lronly)

   drt_voplot,lrs,'/Users/drupke/ifs/esi/plots/makani_esi_ap'+aplab+'_vo_noextcor.eps',$
      errlo=lrsel,errhi=lrseh,scol=['cg1','cg2','cg3']
   drt_voplot,eclrs,'/Users/drupke/ifs/esi/plots/makani_esi_ap'+aplab+'_vo.eps',$
      errlo=eclrsel,errhi=eclrseh,scol=['cg1','cg2','cg3']

   dens_s2 = drt_nebden(10d^lrs['s2'],'s2',err=lrselin['s2'])
   dens_o2 = drt_nebden(10d^lrs['o2'],'o2',err=lrselin['o2'])

   ; fluxes, luminosities, etc.
   ; Planck cosmology
   dist = lumdist(0.459d, H0=67.4d, Omega_M = 0.315d, Lambda0 = 0.685d, /silent)
   haflux = ecfx['Halpha']*1d-16
   hafluxerr = ecfxe['Halpha']*1d-16
   halum = drt_linelum(haflux*1d-3,dist,/ergs,err=hafluxerr*1d-3)
   ;print,halum
   loghalum = alog10(halum[*,0])
   loghalumerr = [[loghalum - alog10(halum[*,0] - halum[*,1])],$
      [alog10(halum[*,0]+halum[*,1])-loghalum]]
   hamass = drt_hiimass(loghalum,dens_s2[*,0],errlum=loghalumerr,$
      errden=dens_s2[*,1:2],/log)
   
   openw,luntxt,'/Users/drupke/ifs/esi/docs/makani_esi_ap'+aplab+'.txt',/get_lun
   
   if iap eq 1 or iap eq 2 or iap eq 5 or iap eq 7 then begin
      ; https://ned.ipac.caltech.edu/level5/Sept12/Calzetti/Calzetti1_2.html
      hasfr = 5.5d-42 * halum
      printf,luntxt,'Halpha SFR (from c1, in Msun/yr): ',hasfr[0,0],'+/-',hasfr[0,1],$
         format='(A0,D0.2,A0,D0.2)'
   endif
   
   printf,luntxt,'Name','ebv','err','logrho_s2','errlo','errhi','logrho_o2','errlo','errhi','hiimass','errlo','errhi',$
      format='(A-10,11A10)
   name = apname+['_1','_2','_T']
   for i=0,n_elements(lrs['ebv'])-1 do begin
      printf,luntxt,name[i],lrs['ebv',i],lrsel['ebv',i],$
         dens_s2[i,0],dens_s2[i,1],dens_s2[i,2],$
         dens_o2[i,0],dens_o2[i,1],dens_o2[i,2],$
         hamass[i,0],hamass[i,1],hamass[i,2],$
         format='(A-10,11D10.4)
   endfor

   printf,luntxt,'Name','z','siginit',format='(A-10,2A10)
   name = apname+['_1','_2']
   for i=0,n_elements(z)-1 do begin
      printf,luntxt,name[i],z[i],sig[i],format='(A-10,D10.6,D10.2)
   endfor

   printf,luntxt,'Name','z','sigma','ebv',format='(A-10,3A10)
   name = apname+['_stel']
   for i=0,n_elements(stel_z)-1 do begin
      printf,luntxt,name[i],stel_z[i],stel_sigma[i],stel_ebv[i],format='(A-10,D10.6,2D10.2)
   endfor

   free_lun,luntxt

;   drt_voplot,lrs,'/Users/drupke/ifs/esi/plots/makani_esi_ap'+aplab+'_vo_noextcor.eps',$
;      errlo=lrsel,errhi=lrseh,scol=['cg1','cg2','cg3']
;   drt_voplot,eclrs,'/Users/drupke/ifs/esi/plots/makani_esi_ap'+aplab+'_vo.eps',$
;      errlo=eclrsel,errhi=eclrseh,scol=['cg1','cg2','cg3']


end
