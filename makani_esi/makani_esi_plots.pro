pro makani_esi_plots

   bad = 1d99
   
   speryr = 24d*3600d*365.25d    ; seconds in a year
   lsun = 3.826d33               ; solar luminosities, erg/s
   msun = 1.989e33               ; solar mass, g
   c_cms = 2.99792d10
   cm_per_pc = 3.0857e18

   plotdir = '/Users/drupke/ifs/esi/plots/'
   plotquantum = 2.5 ; in inches


   ; Planck cosmology
   angdist=1b
   dist = drt_plumdist(0.459d,angdist=angdist)

   kpc_per_as = angdist*1000d/206265d
   esips = 0.1542d
   kpc_per_pix = esips * kpc_per_as

   ; fluxes, line ratios
   restore,'/Users/drupke/specfits/esi/makani/makani.lin.xdr'
   ; redshifts
   restore,'/Users/drupke/specfits/esi/makani/makani.lininit.xdr'
   ; stellar properties
   restore,'/Users/drupke/specfits/esi/makani/makani.cont.xdr'
   ; apertures
   ; skip the first aperture in the fit results because it is the 3.0" spectrum
   apoff = 1
   ; indices of repeated aperture
   irepeat = [6,9]
   readcol,'/Users/drupke/ifs/esi/docs/apertures.csv',apnum,apnames,apradas,$
      appa,apoffpix,aplenpix,s2good1,s2good2,s2goodT,$
      skip=1+apoff,/silent,format='(I,X,A,D,X,D,X,D,D,X,D,D,D)'
   ; calculate aperture radii in kpc
   apradkpc = apradas * kpc_per_as
   ; final aperture is not a separate spectrum, but needs to be overplotted on [OII] image
   naps = n_elements(apnum)-1
   napsp1 = n_elements(apnum)
   ; Is [SII] good for doing density analysis?
   s2good = [[s2good1],[s2good2],[s2goodT]]
   ; multi-component analysis?
   cmult = bytarr(naps)
   cmult[0] = 1b
   ; connect multiple components in plots?
   connect = intarr(naps,3)
   iconnect = 1
   for i=0,naps-1 do begin
      if cmult[i] then begin
         connect[i,*] = iconnect
         iconnect++
      endif
   endfor
   ; symbol sizes
   symsize = dblarr(naps,3) + 1d
   ; Get data
   linelab = 1b
   linelist = ifsf_linelist(/all,/vacuum,linelab=linelab,/quiet)

   ; Compute CVDF velocities
   emlcvdfvel = ifsf_cmpcvdfvals(emlcvdf,emlflx,emlflxerr)

   ; fluxes and errors of lines
   fx = hash()
   fxe = hash()
   ecfx = hash()
   ecfxe = hash()
   ; gas zs and sigmas
   z = dindgen(naps,2) + bad
   sig = dindgen(naps,2) + bad
   ; populate stellar properties
   stel_z = contcube.stel_z[apoff:apoff+naps-1]
   stel_sigma = contcube.stel_sigma[apoff:apoff+naps-1]
   stel_ebv = contcube.stel_ebv[apoff:apoff+naps-1]
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
   ; Hbeta intrinsic SB
   loghblumsb = dindgen(naps,3) + bad
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
   printf,luntxt1,'Ap#','Name','ebv','err','logrho_s2','errlo','errhi','logrho_o2',$
      'errlo','errhi','hiimass','errlo','errhi','SFR','err',$
      format='(A-3,A-15,13A10)
   printf,luntxt2,'Name','z','siginit',format='(A-20,2A10)
   printf,luntxt3,'Name','z','sigma','ebv',format='(A-20,3A10)
   printf,luntxt4,'Name','[OII]/Ha',format='(A-20,1A10)
   
   for iap=0,naps-1 do begin
      iapp1 = iap+apoff
      ; are there multiple lines? Use [OII] to tell
      if emlflx['fc2','[OII]3726+[OII]3729',iapp1] ne bad then ncomp=2 $
      else ncomp=1
      ; populate fluxes and errors of lines
      foreach key,emlflx['fc1'].keys() do begin
         if iap eq 0 then begin
            fx[key] = dindgen(naps,3)+bad
            fxe[key] = dindgen(naps,3)+bad
         endif
         if ncomp eq 2 then begin
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
         ; remove some lines by hand
;         if key eq '[SII]6716+[SII]6731' OR key eq '[SII]6716' or $
;            key eq '[SII]6731' then begin
;               fx[key,5,*] = bad
;               fxe[key,5,*] = bad
;               fx[key,8,*] = bad
;               fxe[key,8,*] = bad
;         endif
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

   ; Set Halpha using Hbeta for aperture 9
   fx['Halpha',8,0:2] = 2.86d*fx['Hbeta',8,0:2]

   ; compute extincted line ratios
   lrs = ifsf_lineratios(fx,fxe,linelist,errlo=lrsel,errhi=lrseh,errlin=lrselin)
  
   for iap=0,naps-1 do begin
      iapp1 = iap+apoff
      ; compute E(B-V)
      foreach key,fx.keys() do begin
         if iap eq 0 then begin
            ecfx[key] = dindgen(naps,3)+bad
            ecfxe[key] = dindgen(naps,3)+bad
         endif
         igd = where(fx[key,iap,*] ne bad AND lrs['ebv',iap,*] ne bad,ctgd)
         if ctgd gt 0 then begin
            ec_tmp = drt_dustcor_ccm(linelist[key],fx[key,iap,igd],$
               lrs['ebv',iap,igd],ebverr=lrselin['ebv',iap,igd],$
               fluxerr=fxe[key,iap,igd])
            ecfx[key,iap,igd] = ec_tmp[*,0]
            ecfxe[key,iap,igd] = ec_tmp[*,1]
         endif
      endforeach
   endfor
   
   ; compute unextincted line ratios
   eclrs = ifsf_lineratios(ecfx,ecfxe,linelist,errlo=eclrsel,errhi=eclrseh,$
      errlin=eclrselin,/lronly)

   igds2 = where(lrs['s2'] ne bad AND s2good[0:naps-1,*])
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
   ; Compute intrinsic Hbeta flux across pixel. Divide by arcsec^2 and by cm^-2/as^-2
   aplenpixarr = rebin(aplenpix[0:naps-1],naps,3)
   loghblumsb[igdha] = alog10(halum) - alog10(aplenpixarr[igdha]*esips) - $
      2d*alog10(kpc_per_as * cm_per_pc * 1d3) - alog10(2.86d)
   ;print,loghblumsb

   ; fill in bad densities with low-den limit
   dens_s2_use = dens_s2[igdha]
   dens_s2errlo_use = dens_s2errlo[igdha]
   dens_s2errhi_use = dens_s2errhi[igdha]
   ibdtmp = where(dens_s2_use eq bad)
   dens_s2_use[ibdtmp] = 0d
   dens_s2errlo_use[ibdtmp] = 0d
   dens_s2errhi_use[ibdtmp] = 0d

   hamass[igdha] = drt_hiimass(loghalum[igdha],dens_s2_use,$
      errlum=[[loghalumerrlo[igdha]],[loghalumerrhi[igdha]]],$
      errden=[[dens_s2errlo_use],[dens_s2errhi_use]],/log,$
      masserr=hamasserr_tmp)
   hamasserrlo[igdha] = hamasserr_tmp[*,0]
   hamasserrhi[igdha] = hamasserr_tmp[*,1]
      
   ; https://ned.ipac.caltech.edu/level5/Sept12/Calzetti/Calzetti1_2.html
   hasfr[igdha] = 7.9d-42 * halum
   hasfrerr[igdha] = 7.9d-42 * halumerr

   for iap=0,naps-1 do begin
      if cmult[iap] then name = apnames[iap]+['_1','_2','_T'] $
      else name = apnames[iap]+['_T']
      for i=0,n_elements(name)-1 do begin
         dens_s2_use = dens_s2[iap,i]
         if dens_s2_use eq bad then dens_s2_use = 0d
         dens_s2errlo_use = dens_s2errlo[iap,i]
         if dens_s2errlo_use eq bad then dens_s2errlo_use = 0d
         dens_s2errhi_use = dens_s2errhi[iap,i]
         if dens_s2errhi_use eq bad then dens_s2errhi_use = 0d
         printf,luntxt1,iap+1,name[i],lrs['ebv',iap,i],lrsel['ebv',iap,i],$
            dens_s2_use,dens_s2errlo_use,dens_s2errhi_use,$
            dens_o2[iap,i],dens_o2errlo[iap,i],dens_o2errhi[iap,i],$
            hamass[iap,i],hamasserrlo[iap,i],hamasserrhi[iap,i],$
            hasfr[iap,i],hasfrerr[iap,i],$
            format='(I-3,A-15,13D10.4)
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


   ; ------------
   ; VO plots
   ; ------------

;   ; use sigma to get symbol size
;   sigsym = [[sig],[emlcvdfvel['vsig','[OII]3726+[OII]3729',apoff:naps-1+apoff]]]
;   igdsig = where(sigsym ne bad)
;   symsize[igdsig] = alog10(sigsym[igdsig]/min(sigsym[igdsig]))+1d
   ; use radius to get symbol size
   radsym = rebin(apradkpc[0:naps-1],naps,3)
   vopsym = intarr(naps,3)+16
   vopsym[*,0] = 9
   vopsym[*,1] = 22
   symsize = ((max(radsym) - radsym)/max(radsym)*1.75) + 0.5d
;   ; use [OII] flux to get symbol size
;   o2sym = fx['[OII]3726+[OII]3729',*,*]
;   igdo2sym = where(o2sym ne bad)
;   maxo2sym = max(o2sym[igdo2sym])
;   symsize = o2sym/maxo2sym*2d + 1d
   ; use radius to get color
   ; http://www.idlcoyote.com/color_tips/brewer.html
   ;cgLoadCT, 0, NCOLORS=16, BOTTOM=0, /BREWER
   ;color = byte((rebin(apradkpc[1:naps],naps,3)/max(apradkpc[1:naps])*12))+4
   tvlct,transpose([228,26,28]),101
   tvlct,transpose([55,126,184]),102
   tvlct,transpose([77,175,74]),103
   tvlct,transpose([152,78,163]),104
   tvlct,transpose([255,127,0]),105
   tvlct,transpose([255,255,51]),106
   tvlct,transpose([166,86,40]),107
   tvlct,transpose([247,129,191]),108
   tvlct,transpose([153,153,153]),109
   apcolor=rebin(bindgen(naps)+101,naps,3)
   ; components we don't want to plot
   mask = bytarr(naps,3)+1b
   for i=0,naps-1 do $
      if not cmult[i] then mask[i,0:1] = 0b
   drt_voplot,lrs,plotdir+'makani_esi_vo_noextcor.eps',$
      errlo=lrsel,errhi=lrseh,connect=connect,symsize=symsize,color=apcolor,$
      mask=mask,psym=vopsym
   drt_voplot,eclrs,plotdir+'makani_esi_vo.eps',$
      errlo=eclrsel,errhi=eclrseh,connect=connect,symsize=symsize,$
      color=apcolor,mask=mask,psym=vopsym

   ; -----------------------------
   ; [OII] map with slits overlaid
   ;-------------------------------
   ;
   ; image geometry / coordinates
   kcwisamplefac = 10
   kcwi_center_nuclei = [30d,35d]
   kcwips = 0.29d
   kcwi_kpc_per_pix =  kcwips * kpc_per_as

   cubedir = '/Users/drupke/ifs/kcwi/cubes/hizeaj2118/'
   ; -300 to +300 km/s, units 10^-16 erg/s/cm^-2
   o2map_core = readfits(cubedir+$
      'hizeaj2118o2_m300_to_p300kms_stelcontsub.fits')
   ; -300 to +300 km/s, units 10^-16 erg/s/cm^-2, subsampled by factor of 10 through interpolation
   o2map_core_ss10 = readfits(cubedir+$
      'hizeaj2118o2_m300_to_p300kms_ss10_stelcontsub.fits')
   ; convert to surface brightness
   o2map_core /= kcwips^2d
   o2map_core_ss10 /= kcwips^2d
   maxo2map_core = max(o2map_core)
   ;print,maxo2map_core
   iclip = where(o2map_core_ss10 le 1.5d-2*maxo2map_core,ctclip)
   if ctclip gt 0 then o2map_core_ss10[iclip] = 0d
   
   ; more image geometry / coordinates
   size_tmp = size(o2map_core)
   kcwidx = size_tmp[1]
   kcwidy = size_tmp[2]
   kcwi_center_axes = [double(kcwidx)/2d,double(kcwidy)/2d]+0.5d
   kcwi_xran_kpc = double([-(kcwi_center_axes[0]-0.5),kcwidx-(kcwi_center_axes[0]-0.5)]) $
      * kcwi_kpc_per_pix
   kcwi_yran_kpc = double([-(kcwi_center_axes[1]-0.5),kcwidy-(kcwi_center_axes[1]-0.5)]) $
      * kcwi_kpc_per_pix
   kcwi_center_nuclei_kpc_x = (kcwi_center_nuclei[0,*] - kcwi_center_axes[0]) $
      * kcwi_kpc_per_pix
   kcwi_center_nuclei_kpc_y = (kcwi_center_nuclei[1,*] - kcwi_center_axes[1]) $
      * kcwi_kpc_per_pix

   
   cgps_open,plotdir+'makani-o2-esi-slits.eps',charsize=1.5,$
      /encap,/inches,xs=plotquantum*2.5,ys=plotquantum*2.25,/qui,/nomatch
   pos = cglayout([1,1],oxmar=[1,6],oymar=[5,1],$
      xgap=0,ygap=0,aspect=double(kcwidy)/double(kcwidx),unit=!D.X_PX_CM/3.0)

   zran = [0.015,1]*maxo2map_core
   ncbdiv = 4
   dzran = zran[1]-zran[0]
   cbform = '(D0.3)'
   cbpos=[pos[2,0],pos[1,0],pos[2,0]+0.02,pos[3,0]]

   mapscl = cgimgscl(o2map_core_ss10,minval=zran[0],maxval=zran[1],$
      stretch=4,constant=120d)
;   mapscl = cgimgscl(o2map_core,minval=zran[0],maxval=zran[1],$
;      stretch=4,constant=60d)
   ;loadcv,121,/noqual,/silent
   cgloadct,0
   cgimage,mapscl,pos=pos[*,0],/keep,opos=truepos
   ;  line contours
   cgplot,[0],xsty=5,ysty=5,pos=truepos,/nodata,/noerase,$
      xran=[0,kcwidx*kcwisamplefac],yran=[0,kcwidy*kcwisamplefac]
;   cgplot,[0],xsty=5,ysty=5,pos=truepos,/nodata,/noerase,$
;      xran=[0,kcwidx],yran=[0,kcwidy]
   contourlevels = [0.02,0.04,0.08,0.16,0.32,0.64]*2.54
   cgcontour,o2map_core_ss10,dindgen(kcwidx*kcwisamplefac)+0.5,$
      dindgen(kcwidy*kcwisamplefac)+0.5,$
      /overplot,c_thick=1,c_color='white',$
      levels=contourlevels,label=0
;   cgcontour,o2map_core,dindgen(kcwidx)+0.5,$
;      dindgen(kcwidy)+0.5,$
;      /overplot,c_thick=1,c_color='white',$
;      levels=contourlevels,label=0
   cgplot,[0],xsty=5,ysty=5,pos=truepos,/nodata,/noerase,$
      xran=[0,kcwidx],yran=[0,kcwidy]
   ; slits
   ; PA, East of North
   slit_pa = [45d,0d]
   ; Center of slit in KCWI pixel coordinates
   slit_pos = [[kcwi_center_nuclei-0.5d],$
      [kcwi_center_nuclei-0.5d + [3.3d/kcwips,0]]]
   ; slit dimensions in KCWI pixels
   slitlength = 20d / kcwips
   slitwidth = 1d / kcwips
   halflength = slitlength/2d
   halfwidth = slitwidth/2d
   ; slit coordinates in KCWI pixels, in coordinates where [0,0] is center of slit
   ; 4 corners + repeat of first corner
   slitcoord=[[-halfwidth,-halflength],[-halfwidth,halflength],$
      [halfwidth,halflength],[halfwidth,-halflength],[-halfwidth,-halflength]]
   ; loop over slits
   for i=0,n_elements(slit_pa)-1 do begin
      sinang = sin(slit_pa[i]*!DPi/180d)
      cosang = cos(slit_pa[i]*!DPi/180d)
      ; rotate coordinates counterclockwise through angle = PA
      rotmat = [[cosang,sinang],[-sinang,cosang]]
      slitcoordrot = slitcoord
      for j=0,n_elements(slitcoord[0,*])-1 do $
         slitcoordrot[*,j] = rotmat # slitcoord[*,j]
      ; ... and shift
      slitcoordrotshift = rebin(slit_pos[*,i],2,5) + slitcoordrot
      cgpolygon,slitcoordrotshift,color='White'
   endfor
   ; apertures
   ; PA, East of North
   ; Location of slit center in KCWI pixel coordinates
   appos45 = kcwi_center_nuclei-0.5d
   appos0 = kcwi_center_nuclei-0.5d + [3.3d/kcwips,0]
   appos = []
   for i=0,n_elements(appa)-1 do begin
      if appa[i] eq 45 then appos = [[appos],[appos45]] $
      else appos = [[appos],[appos0]]
   endfor
   ; center of aperture in unrotated KCWI pixels
   apcent = [[transpose(dblarr(napsp1))],[transpose(apoffpix)]*esips/kcwips]
   ; loop over aps
   ;cgLoadCT, 0, NCOLORS=16, BOTTOM=0, /BREWER
   ;color = byte((rebin(apradkpc[1:nappos],nappos,3)/max(apradkpc[1:nappos])*12))+4
   tvlct,transpose([228,26,28]),101
   tvlct,transpose([55,126,184]),102
   tvlct,transpose([77,175,74]),103
   tvlct,transpose([152,78,163]),104
   tvlct,transpose([255,127,0]),105
   tvlct,transpose([255,255,51]),106
   tvlct,transpose([166,86,40]),107
   tvlct,transpose([247,129,191]),108
   tvlct,transpose([153,153,153]),109
   apcolor=bindgen(napsp1)+101
   ; repeated aperture
   apcolor[irepeat[1]] = apcolor[irepeat[0]]
   ; array to hold galactocentric radii of aperture corners
   apcoord_radius_kpc = dblarr(napsp1,5)
   for i=0,n_elements(appa)-1 do begin
      ; ap dimensions in KCWI pixels
      halfwidth = 1d / kcwips/2d
      halflength = aplenpix[i]*esips / kcwips / 2d
      ; coordinates of aperture corners centered on [0,0]
      apcoord = [[-halfwidth,-halflength],[-halfwidth,halflength],$
         [halfwidth,halflength],[halfwidth,-halflength],[-halfwidth,-halflength]]
      ; coordinates of aperture corners shifted by aperture center
      apcoord = apcoord + rebin(apcent[*,i],2,5)
      sinang = sin(appa[i]*!DPi/180d)
      cosang = cos(appa[i]*!DPi/180d)
      ; rotate coordinates counterclockwise through angle = PA about [0,0]
      rotmat = [[cosang,sinang],[-sinang,cosang]]
      apcoordrot = apcoord
      for j=0,n_elements(apcoord[0,*])-1 do $
         apcoordrot[*,j] = rotmat # apcoord[*,j]
      ; ... and shift by nucleus location
      apcoordrotshift = rebin(appos[*,i],2,5) + apcoordrot
      cgpolygon,apcoordrotshift,color=apcolor[i],thick=8
      ; cgpolygon changes apcoordrotshift
      apcoordrotshift = rebin(appos[*,i],2,5) + apcoordrot
      apcoord_xy_pix = apcoordrotshift - rebin(kcwi_center_nuclei-0.5d,2,5)
      apcoord_radius_pix = sqrt(apcoord_xy_pix[0,*]^2d + apcoord_xy_pix[1,*]^2d)
      apcoord_radius_kpc[i,*] = apcoord_radius_pix*kcwi_kpc_per_pix
   endfor

   ;
   ifsf_plotaxesnuc,kcwi_xran_kpc,kcwi_yran_kpc,kcwi_center_nuclei_kpc_x,$
      kcwi_center_nuclei_kpc_y ;,colorax='white'
   al_legend,string(indgen(naps)+1,format='(I0)'),textcolors=apcolor[0:naps-1],/bottom,/left
;   including aperture 10 (bottom of 0 degree slit) as separate aperture label
;   al_legend,string(indgen(napsp1)+1,format='(I0)'),textcolors=apcolor,/bottom,/left
;   ;ifsf_plotcompass,xarr_kpc,yarr_kpc,carr='Black',hsize=150d,hthick=2d
   ticknames = strarr(ncbdiv+1)
   ticknames[0] = zran[0]
   ticknames[ncbdiv] = zran[1]
   dticks = 256d/ncbdiv
   for i=1,ncbdiv-1 do begin
      ticks_byte = byte(i*dticks)
      ibyte = where(mapscl eq ticks_byte)
;      ticknames[i] = mean(o2map_core[ibyte])
      ticknames[i] = mean(o2map_core_ss10[ibyte])
   endfor
   ticknames = string(ticknames,format=cbform)
   cgloadct,0
   cgcolorBar,position=cbpos,ticknames=ticknames,divisions=ncbdiv,/ver,/right,$
      charsize=1
   cgtext,0.96,0.57,'F ([OII], -300 to +300 km s$\up-1$, 10$\up-16$ erg s$\up-1$ cm$\up-2$ arcsec$\up-2$)',$
      orient=270,/normal,align=0.5,chars=1
   cgtext,0.45,0.02,'$\Delta$x (kpc)',align=0.5,/normal
   cgtext,0.05,0.55,'$\Delta$y (kpc)',align=0.5,/normal,orient=90

   cgps_close


   ; --------------------------
   ;  [OII], Halpha vs. radius
   ; --------------------------

   tvlct,transpose([228,26,28]),101
   tvlct,transpose([55,126,184]),102
   tvlct,transpose([77,175,74]),103
   tvlct,transpose([152,78,163]),104
   tvlct,transpose([255,127,0]),105
   tvlct,transpose([255,255,51]),106
   tvlct,transpose([166,86,40]),107
   tvlct,transpose([247,129,191]),108
   tvlct,transpose([153,153,153]),109

   ; Compute line ratios
   lrlist=hash()
   lrlist['o2ha']=['[OII]3726+[OII]3729','Halpha']
   ; extincted line ratio
   o2ha_ext = ifsf_lineratios(fx,fxe,linelist,lrlist=lrlist,$
      errlo=o2ha_ext_el,errhi=o2ha_ext_eh,errlin=o2ha_ext_elin,/lronly)
   ; unextincted line ratio
   o2ha = ifsf_lineratios(ecfx,ecfxe,linelist,lrlist=lrlist,$
      errlo=o2ha_el,errhi=o2ha_eh,errlin=o2ha_elin,/lronly)
   ; [OII] extincted / Halpha intrinsic
   o2ha_exto2 = o2ha_ext[*]
   for iap=0,naps-1 do begin      
      igd = where(fx['Halpha',iap,*] ne bad AND lrs['ebv',iap,*] ne bad,ctgd)
      if ctgd gt 0 then $
         o2ha_exto2['o2ha',iap,igd] += alog10(fx['Halpha',iap,igd]) $
            - alog10(ecfx['Halpha',iap,igd])
   endfor

   ; plots

   cgps_open,plotdir+'makani_esi_o2ha_v_rad.eps',$
      xsize=6,ysize=9,/inches,/encap,charsize=1,default_thick=1,/qui,/nomatch
   pos = cglayout([1,4],oxmar=[8,1],oymar=[5,1],xgap=0,ygap=0)

   ; to overplot models
   xmod = dindgen(100d)/99d*35d

   ; extincted fluxes 
   x = apradkpc[0:naps-1]
   xel = x-min(apcoord_radius_kpc,dim=2)
   apradkpc_lo = xel
   xel = xel[0:naps-1]
   xeh = max(apcoord_radius_kpc,dim=2)-x
   apradkpc_hi = xeh
   xeh = xeh[0:naps-1]
   xran = [-5,40]
   xtit = 'Radius (kpc)'
   key1 = '[OII]3726+[OII]3729'
   ; Normalize to surface brightness
   y = alog10(fx[key1,*,2]) - alog10(aplenpix*esips)
   yel = alog10(fx[key1,*,2]) - alog10(fx[key1,*,2] - fxe[key1,*,2])
   yeh = alog10(fx[key1,*,2] + fxe[key1,*,2]) - alog10(fx[key1,*,2])
   key2 = 'Halpha'
   ; Normalize to surface brightness
   yh = alog10(fx[key2,*,2]) - alog10(aplenpix*esips)
   yhel = alog10(fx[key2,*,2]) - alog10(fx[key2,*,2] - fxe[key2,*,2])
   yheh = alog10(fx[key2,*,2] + fxe[key2,*,2]) - alog10(fx[key2,*,2])
   yran = [-2.75,1.5]
   ytit = 'log[$\Sigma$$\downobs$/10$\up-16$ erg s$\up-1$ cm$\up-2$ as$\up-2$]'
   xoff = (randomu(seed,naps*3)-0.5d)*1.5d
   ;xoff = dblarr(naps)
   xoff[0] = 0d

   readcol,'/Users/drupke/ifs/kcwi/maps/hizeaj2118/rb1/hizeaj2118.o2radprof.txt',$
      o2radprof_rad,o2radprof_sb,o2radprof_el,o2radprof_eh,skip=1,/silent,$
      format='(D,D,D,D)'

   cgplot,[0],xran=xran,yran=yran,ytit=ytit,xtickf='(A1)',$
      /xsty,/ysty,/nodata,pos=pos[*,0]
   cgoplot,o2radprof_rad,o2radprof_sb+16d,err_ylo=o2radprof_el,$
      err_yhi=o2radprof_eh,/err_clip,err_width=0d
   for i=0,n_elements(x)-1 do begin
      cgoplot,x[i]+xoff[i],y[i],psym=6,symsize=symsize[i],color='Blue',$
         err_xlo=xel[i],err_xhi=xeh[i],err_ylo=yel[i],err_yhi=yeh[i],$
         err_color='Grey',err_thick=2,/err_clip,err_width=0d
      cgoplot,x[i]+xoff[i],yh[i],psym=16,symsize=symsize[i],color='Red',$
         err_xlo=xel[i],err_xhi=xeh[i],err_ylo=yhel[i],err_yhi=yheh[i],$
         err_color='Grey',err_thick=2,/err_clip,err_width=0d
      cgoplot,[x[i]+xoff[i],x[i]+xoff[i]],[y[i],yh[i]]
   endfor
   al_legend,['[OII]','H$\alpha$'],color=['Blue','Red'],psym=[6,16],/top,/right

   ymod = xmod^2d

   ; ext vs rad
   yran = [-0.1,1.25]
   y = lrs['ebv',*,2]
   yel = lrsel['ebv',*,2]
   yeh = lrseh['ebv',*,2]
   ytit = 'E(B-V)'
   cgplot,[0],xran=xran,yran=yran,ytit=ytit,xtickf='(A1)',$
      /xsty,/ysty,/nodata,pos=pos[*,1],/noerase
   for i=0,n_elements(x)-1 do $
      cgoplot,x[i]+xoff[i],y[i],psym=16,symsize=symsize[i],color=apcolor[i],$
      err_xlo=xel[i],err_xhi=xeh[i],err_ylo=yel[i],err_yhi=yeh[i],$
      err_color='Grey',err_thick=2,/err_clip,err_width=0d

   yran = [-0.75,1.25]
   ye = o2ha_ext['o2ha',*,2]
   yeel = o2ha_ext_el['o2ha',*,2]
   yeeh = o2ha_ext_eh['o2ha',*,2]
   y = o2ha['o2ha',*,2]
   yel = o2ha_el['o2ha',*,2]
   yeh = o2ha_eh['o2ha',*,2]
   ytit = 'log([OII]/H$\alpha$)
   cgplot,[0],xran=xran,yran=yran,ytit=ytit,xtickf='(A1)',$
      /xsty,/ysty,/nodata,pos=pos[*,2],/noerase
   for i=0,n_elements(x)-1 do begin
      if i eq n_elements(x)-1 then begin
         cgoplot,x[i]+xoff[i],y[i],psym=5,symsize=symsize[i]*1.5,color='Blue',$
            err_xlo=xel[i],err_xhi=xeh[i],err_yhi=yeeh[i],$
            err_color='Grey',err_thick=2,/err_clip,err_width=0d
         cgoplot,x[i]+xoff[i],ye[i],psym=16,symsize=symsize[i],color='Red',$
            err_xlo=xel[i],err_xhi=xeh[i],err_ylo=yeel[i],err_yei=yeeh[i],$
            err_color='Grey',err_thick=2,/err_clip,err_width=0d
      endif else begin
         cgoplot,x[i]+xoff[i],y[i],psym=6,symsize=symsize[i],color='Blue',$
            err_xlo=xel[i],err_xhi=xeh[i],err_ylo=yel[i],err_yhi=yeh[i],$
            err_color='Grey',err_thick=2,/err_clip,err_width=0d
         cgoplot,x[i]+xoff[i],ye[i],psym=16,symsize=symsize[i],color='Red',$
            err_xlo=xel[i],err_xhi=xeh[i],err_ylo=yeel[i],err_yei=yeeh[i],$
            err_color='Grey',err_thick=2,/err_clip,err_width=0d
      endelse
      cgoplot,[x[i]+xoff[i],x[i]+xoff[i]],[y[i],ye[i]]
   endfor
   al_legend,['intrinsic','observed'],color=['Blue','Red'],psym=[6,16],/bottom,/right

   yran = [-1.5,1.5]
   y = o2ha_exto2['o2ha',*,2]
   yel = o2ha_ext_el['o2ha',*,2]
   yeh = o2ha_ext_eh['o2ha',*,2]
   ytit = 'log( [OII], observed / H$\alpha$, intrinsic )'
   ; Fit [OII]obs/Haint
   xfit = x[0:naps-1]
   yfit = y[0:naps-1]
   nfit = naps-1
   yfiterr = (yel[0:naps-1]+yeh[0:naps-1])/2d
   fitpower = 1d
   ;o2ha_exto2_pars =  mpfitexpr('P[0]+P[1]*X',xfit,yfit,yfiterr,/quiet)
   o2ha_exto2_pars =  mpfitexpr('P[0]+P[1]*X',xfit,yfit,yfiterr,/quiet)
   ;xfitlo = min(xfit)
   ;xfithi = max(xfit)
   ;ixfitlo = where(xfit eq xfitlo)
   ;ixfithi = where(xfit eq xfithi)
   ;yfitlo = o2ha_exto2_pars[0]+xfit[ixfitlo[0]]*o2ha_exto2_pars[1]
   ;yfithi = o2ha_exto2_pars[0]+xfit[ixfithi[0]]*o2ha_exto2_pars[1]
   ;yfitlo = o2ha_exto2_pars[0]+$ ;xfit[ixfitlo[0]]*o2ha_exto2_pars[1] + $
   ;   xfit[ixfithi[0]]^2d*o2ha_exto2_pars[1]
   ;yfithi = o2ha_exto2_pars[0]+$ ;xfit[ixfithi[0]]*o2ha_exto2_pars[1]+ $
   ;   xfit[ixfithi[0]]^2d*o2ha_exto2_pars[1]
   yfitmod = o2ha_exto2_pars[0]+xfit^fitpower*o2ha_exto2_pars[1]

   cgplot,[0],xran=xran,yran=yran,ytit=ytit,xtit=xtit,$
      /xsty,/ysty,/nodata,pos=pos[*,3],/noerase
   ;cgoplot,xfit,o2ha_exto2_pars[0]+xfit*o2ha_exto2_pars[1],thick=2
   cgoplot,xmod,o2ha_exto2_pars[0]+ $ ;xmod*o2ha_exto2_pars[1]+$
      xmod^fitpower*o2ha_exto2_pars[1],thick=2

   ; RMS for case of fit to every point
   ;o2ha_exto2_rms = sqrt(mean((yfit-(o2ha_exto2_pars[0]+xfit*o2ha_exto2_pars[1]))^2d))
   ;cgoplot,xfit,o2ha_exto2_pars[0]+xfit*o2ha_exto2_pars[1]+o2ha_exto2_rms,linesty=1
   ;cgoplot,xfit,o2ha_exto2_pars[0]+xfit*o2ha_exto2_pars[1]-o2ha_exto2_rms,linesty=1
   o2ha_exto2_rms = sqrt(mean((yfit-yfitmod)^2d))
   cgoplot,xmod,o2ha_exto2_pars[0]+xmod^fitpower*o2ha_exto2_pars[1]+o2ha_exto2_rms,linesty=1
   cgoplot,xmod,o2ha_exto2_pars[0]+xmod^fitpower*o2ha_exto2_pars[1]-o2ha_exto2_rms,linesty=1

   ;cgoplot,[0,xfitlo],[yfitlo,yfitlo],thick=2
   ;cgoplot,[xfithi,xran[1]],[yfithi,yfithi],thick=2
   for i=0,n_elements(x)-1 do $
      cgoplot,x[i]+xoff[i],y[i],psym=16,symsize=symsize[i],color=apcolor[i],$
         err_xlo=xel[i],err_xhi=xeh[i],err_ylo=yel[i],err_yhi=yeh[i],$
         err_color='Grey',err_thick=2,/err_clip,err_width=0d
;   cgtext,-3d,1.0,'log([OII],obs/Ha,int) = '+$
;      string(o2ha_exto2_pars[0],format='(D0.2)')+' + '+$
;      string(o2ha_exto2_pars[1],format='(D0.3)')+' (r/kpc)'
   cgtext,-3d,1.0,'log([OII],obs/Ha,int) = '+$
      string(o2ha_exto2_pars[0],format='(D0.2)')+' + '+$
      string(o2ha_exto2_pars[1],format='(E0.2)')+' (r/kpc)' ;$\up2$'
   cgtext,-3d,0.7,'RMS = '+string(o2ha_exto2_rms,format='(D0.2)')+' dex'
 

   cgps_close

;
;   ; --------------------------
;   ;  [OII], Halpha vs. [OII] flux
;   ; --------------------------
;
;   ; plots
;
;   cgps_open,plotdir+'makani_esi_o2ha_v_o2.eps',$
;      xsize=6,ysize=9,/inches,/encap,charsize=1,default_thick=1,/qui,/nomatch
;   pos = cglayout([1,4],oxmar=[8,1],oymar=[5,1],xgap=0,ygap=0)
;
;   ; extincted fluxes
;   xran = [-1.75,1.]
;   xtit = 'log[$\Sigma$$\downobs$([OII])/10$\up-16$ erg s$\up-1$ cm$\up-2$ as$\up-2$]'
;   key1 = '[OII]3726+[OII]3729'
;   y = alog10(fx[key1,*,2]) - alog10(aplenpix*esips)
;   yel = y - alog10(fx[key1,*,2] - fxe[key1,*,2])
;   yeh = alog10(fx[key1,*,2] + fxe[key1,*,2]) - y
;   x = y
;   xel = yel
;   xeh = yeh
;   key2 = 'Halpha'
;   yh = alog10(fx[key2,*,2]) - alog10(aplenpix*esips)
;   yhel = yh - alog10(fx[key2,*,2] - fxe[key2,*,2])
;   yheh = alog10(fx[key2,*,2] + fxe[key2,*,2]) - yh
;   yran = [-1.75,1.5]
;   ytit = 'log[$\Sigma$$\downobs$/10$\up-16$ erg s$\up-1$ cm$\up-2$ as$\up-2$]'
;   cgplot,[0],xran=xran,yran=yran,ytit=ytit,xtickf='(A1)',$
;      /xsty,/ysty,/nodata,pos=pos[*,0]
;   for i=0,n_elements(x)-1 do begin
;      cgoplot,x[i],y[i],psym=16,symsize=symsize[i],color='Blue',$
;         err_xlo=xel[i],err_xhi=xeh[i],err_ylo=yel[i],err_yhi=yeh[i],$
;         err_color='Grey',err_thick=2,/err_clip,err_width=0d
;      cgoplot,x[i],yh[i],psym=16,symsize=symsize[i],color='Red',$
;         err_xlo=xel[i],err_xhi=xeh[i],err_ylo=yhel[i],err_yhi=yheh[i],$
;         err_color='Grey',err_thick=2,/err_clip,err_width=0d
;      cgoplot,[x[i],x[i]],[y[i],yh[i]]
;   endfor
;   al_legend,['[OII]','H$\alpha$'],color=['Blue','Red'],psym=16,/bottom,/right
;
;   ; ext
;   yran = [-0.25,1.25]
;   y = lrs['ebv',*,2]
;   yel = lrsel['ebv',*,2]
;   yeh = lrseh['ebv',*,2]
;   ytit = 'E(B-V)'
;   cgplot,[0],xran=xran,yran=yran,ytit=ytit,xtickf='(A1)',$
;      /xsty,/ysty,/nodata,pos=pos[*,1],/noerase
;   for i=0,n_elements(x)-1 do $
;      cgoplot,x[i],y[i],psym=16,symsize=symsize[i],color='Blue',$
;      err_xlo=xel[i],err_xhi=xeh[i],err_ylo=yel[i],err_yhi=yeh[i],$
;      err_color='Grey',err_thick=2,/err_clip,err_width=0d
;
;   yran = [-0.75,1]
;   ye = o2ha_ext['o2ha',*,2]
;   yeel = o2ha_ext_el['o2ha',*,2]
;   yeeh = o2ha_ext_eh['o2ha',*,2]
;   y = o2ha['o2ha',*,2]
;   yel = o2ha_el['o2ha',*,2]
;   yeh = o2ha_eh['o2ha',*,2]
;   ytit = 'log([OII]/H$\alpha$)
;   cgplot,[0],xran=xran,yran=yran,ytit=ytit,xtickf='(A1)',$
;      /xsty,/ysty,/nodata,pos=pos[*,2],/noerase
;   for i=0,n_elements(x)-1 do begin
;      cgoplot,x[i],y[i],psym=16,symsize=symsize[i],color='Blue',$
;         err_xlo=xel[i],err_xhi=xeh[i],err_ylo=yel[i],err_yhi=yeh[i],$
;         err_color='Grey',err_thick=2,/err_clip,err_width=0d
;      cgoplot,x[i],ye[i],psym=16,symsize=symsize[i],color='Red',$
;         err_xlo=xel[i],err_xhi=xeh[i],err_ylo=yeel[i],err_yei=yeeh[i],$
;         err_color='Grey',err_thick=2,/err_clip,err_width=0d
;      cgoplot,[x[i],x[i]],[y[i],ye[i]]
;   endfor
;   al_legend,['intrinsic','observed'],color=['Blue','Red'],psym=16,/top,/right
;
;   yran = [-1.5,0.75]
;   y = o2ha_exto2['o2ha',*,2]
;   yel = o2ha_el['o2ha',*,2]
;   yeh = o2ha_eh['o2ha',*,2]
;   ytit = 'log( [OII], observed / H$\alpha$, intrinsic )
;   cgplot,[0],xran=xran,yran=yran,ytit=ytit,xtit=xtit,$
;      /xsty,/ysty,/nodata,pos=pos[*,3],/noerase
;   for i=0,n_elements(x)-1 do $
;      cgoplot,x[i],y[i],psym=16,symsize=symsize[i],color='Blue',$
;         err_xlo=xel[i],err_xhi=xeh[i],err_ylo=yel[i],err_yhi=yeh[i],$
;         err_color='Grey',err_thick=2,/err_clip,err_width=0d
;
;   cgps_close
;



   ; --------------------------
   ;  Density vs. radius
   ; --------------------------

   cgps_open,plotdir+'makani_esi_s2ne_v_rad.eps',$
      xsize=7.5,ysize=7.5,/inches,/encap,charsize=1,default_thick=1,/qui,/nomatch
   x = rebin(apradkpc[0:naps-1],naps,3)
   xel = x - rebin((min(apcoord_radius_kpc,dim=2))[0:naps-1],naps,3)
   xeh = rebin((max(apcoord_radius_kpc,dim=2))[0:naps-1],naps,3) - x
   xran = [-5,30]
   xtit = 'Radius (kpc)'
   y = dens_s2
   yel = dens_s2errlo
   yeh = dens_s2errhi
   yran = [-1,6]
   ytit = 'log(n$\downe$/cm$\up-3$) ([SII] ratio)'
   xoff = (randomu(seed,naps*3)-0.5d)*2d
   apcolor=rebin(bindgen(naps)+101,naps,3)
   cgplot,[0],xran=xran,yran=yran,ytit=ytit,xtit=xtit,$
      /xsty,/ysty,/nodata
   for i=0,n_elements(x)-1 do begin
      if y[i] ne bad then $
         cgoplot,x[i]+xoff[i],y[i],psym=16,symsize=symsize[i],color=apcolor[i],$
         err_xlo=xel[i],err_xhi=xeh[i],err_ylo=yel[i],err_yhi=yeh[i],$
         err_color='Grey',err_thick=2,/err_clip,err_width=0d
   endfor
   if keyword_set(connect) then begin
      for j=1,max(connect) do begin
         indj = where(connect eq j)
         cgoplot,x[indj]+xoff[indj],y[indj]
      endfor
   endif
   cgps_close

   ; --------------------------
   ;  Neutral tracers vs. radius
   ; --------------------------

   tvlct,transpose([27,158,119]),101
   tvlct,transpose([217,95,2]),102
   tvlct,transpose([117,112,179]),103
   ;   tvlct,transpose([152,78,163]),104
   ;   tvlct,transpose([255,127,0]),105
   ;   tvlct,transpose([255,255,51]),106
   ;   tvlct,transpose([166,86,40]),107
   ;   tvlct,transpose([247,129,191]),108
   ;   tvlct,transpose([153,153,153]),109

   ; plots

   cgps_open,plotdir+'makani_esi_neutral_v_rad.eps',$
      xsize=7.5,ysize=7.5,/inches,/encap,charsize=1,default_thick=1,/qui,/nomatch
   ;xran = [-0.1,2.4]
   xran = [-5,60]
   yran = [-2,1]
   ;xtit = 'log(radius / kpc)'
   xtit = 'Radius (kpc)'
   ytit = 'log[$\Sigma$$\downobs$/10$\up-16$ erg s$\up-1$ cm$\up-2$ as$\up-2$]'

   readcol,'/Users/drupke/ifs/kcwi/maps/hizeaj2118/rb1/hizeaj2118.o2radprof.txt',$
      o2radprof_rad,o2radprof_sb,o2radprof_el,o2radprof_eh,skip=1,/silent,$
      format='(D,D,D,D)'
   readcol,'/Users/drupke/ifs/kcwi/maps/hizeaj2118/rb1/hizeaj2118.mg2radprof.txt',$
      mg2radprof_rad,mg2radprof_sb,mg2radprof_el,mg2radprof_eh,skip=1,/silent,$
      format='(D,D,D,D)'
   nadflux = [0.32d,0.20d]
   nadfluxerr = [0.02d,0.02d]
   nadsb = alog10(nadflux) - alog10(aplenpix[0:1]*esips)
   nadsblo = alog10(nadflux) - alog10(nadflux - nadfluxerr)
   nadsbhi = alog10(nadflux + nadfluxerr) - alog10(nadflux)
   ;xlog = [0d,alog10(x[1])]
   ;xloglo = [0d,alog10(x[1])-alog10(x[1]-xel[1])]
   ;xloghi = [alog10(xeh[0]),alog10(x[1]+xeh[1])-alog10(x[1])]

   cgplot,[0],xran=xran,yran=yran,xtit=xtit,ytit=ytit,/xsty,/ysty,/nodata
   ;   cgoplot,alog10(o2radprof_rad),o2radprof_sb+16d,err_ylo=o2radprof_el,$
   ;      err_yhi=o2radprof_eh,/err_clip,err_width=0d,color='Blue'
   ;   cgoplot,alog10(o2radprof_rad),o2radprof_sb+16d - (o2ha_exto2_pars[0]+ $
   ;      o2radprof_rad^fitpower*o2ha_exto2_pars[1]),err_ylo=0.3d,err_yhi=0.3d,$
   ;      /err_clip,err_width=0d,thick=2,color='Magenta'
   ;   cgoplot,alog10(mg2radprof_rad),mg2radprof_sb+16d,err_ylo=mg2radprof_el,$
   ;      err_yhi=mg2radprof_eh,/err_clip,err_width=0d,color='Red'
   ;   cgoplot,xlog,nadsb,psym=6,symsize=2,color='Green',$
   ;      err_xlo=xloglog,err_xhi=xloghi,err_ylo=nadsblo,err_yhi=nadsbhi,$
   ;      err_color='Grey',err_thick=2,/err_clip,err_width=0d
   ;   cgoplot,[alog10(9d),alog10(9d)],yran,linesty=1,color='Red'
   ;   cgoplot,[alog10(17d),alog10(17d)],yran,linesty=1,color='Blue'
   ;   cgoplot,[alog10(25d),alog10(25d)],yran,linesty=1,color='Green'
   ;   cgoplot,[alog10(50d),alog10(50d)],yran,linesty=1,color='Black'
   cgoplot,o2radprof_rad,o2radprof_sb+16d,err_ylo=o2radprof_el,$
      err_yhi=o2radprof_eh,/err_clip,err_width=0d,color=101
   cgoplot,[17d,17d],[-0.3,yran[1]],linesty=2,color=101
   cgtext,-4,0.6,'[OII]',color=101
   cgtext,16,0,'[OII] 1/2-light rad',orient=90,color=101
   ;   cgoplot,o2radprof_rad,o2radprof_sb+16d - (o2ha_exto2_pars[0]+ $
   ;      o2radprof_rad^fitpower*o2ha_exto2_pars[1]),err_ylo=0.3d,err_yhi=0.3d,$
   ;      /err_clip,err_width=0d,thick=2,color='Magenta'
   cgoplot,mg2radprof_rad,mg2radprof_sb+16d,err_ylo=mg2radprof_el,$
      err_yhi=mg2radprof_eh,/err_clip,err_width=0d,color=103
   cgoplot,[9d,9d],[yran[0],-1.1],linesty=2,color=103
   cgtext,8,-1.8,'MgII 1/2-light rad',orient=90,color=103
   cgtext,-4,0.1,'MgII',color=103
   cgoplot,x[0:1],nadsb,psym=6,symsize=2,color=102,$
      err_xlo=xel[0:1],err_xhi=xeh[0:1],err_ylo=nadsblo,err_yhi=nadsbhi,$
      err_color='Grey',err_thick=2,/err_clip,err_width=0d
   cgtext,-4,-0.6,'NaI',color=102
   cgoplot,[25d,25d],yran,linesty=1,color=102
   cgarrow,25,0,35,0,/data,thick=4,color=102,hthick=2
   cgtext,27,0.2,'No dust or neutral/molecular',color=102
   cgtext,28,0.1,'gas detected at r > 20-25 kpc',color=102
   cgoplot,[50d,50d],[yran[0],-1.5],color='Black'
   cgtext,44,-1.8,'Wind'
   cgtext,44,-1.9,'Radius'
   cgoplot,xran,[-1.6,-1.6],linesty=1
   cgtext,30,-1.7,'detection limit'

   ;   cgplot,[0],xsty=5,ysty=5,$
   ;      /nodata,/noerase,xran=xran,yran=yran+1.5d
   ;   x = rebin(apradkpc[0:naps-1],naps,3)
   ;   xel = x - rebin((min(apcoord_radius_kpc,dim=2))[0:naps-1],naps,3)
   ;   xeh = rebin((max(apcoord_radius_kpc,dim=2))[0:naps-1],naps,3) - x
   ;   y = dens_s2
   ;   yel = dens_s2errlo
   ;   yeh = dens_s2errhi
   ;   ;ytit = 'log(n$\downe$/cm$\up-3$) ([SII] ratio)'
   ;   xoff = rebin((randomu(seed,naps*3)-0.5d)*2d,naps,3)
   ;   xoff[0,*] = 0d
   ;   apcolor=rebin(bindgen(naps)+101,naps,3)
   ;   for i=0,n_elements(x)-1 do begin
   ;      if y[i] ne bad then $
   ;         cgoplot,x[i]+xoff[i],y[i],psym=16,symsize=symsize[i],color=apcolor[i],$
   ;            err_xlo=xel[i],err_xhi=xeh[i],err_ylo=yel[i],err_yhi=yeh[i],$
   ;            err_color='Grey',err_thick=2,/err_clip,err_width=0d
   ;   endfor
   ;   if keyword_set(connect) then begin
   ;      for j=1,max(connect) do begin
   ;         indj = where(connect eq j)
   ;         cgoplot,x[indj]+xoff[indj],y[indj]
   ;      endfor
   ;   endif


   cgps_close
   ;
   ;


   ;------------
   ; Spectra
   ;------------
   linoth = strarr(2,6)
   linoth[0,0] = ['[OII]3726','']
   linoth[*,1] = ['[NeIII]3967','Hepsilon']
   linoth[*,2] = ['[OIII]4959','[OIII]5007']
   linoth[*,3] = ['[OI]6364','']
   linoth[*,4] = ['Halpha','[NII]6583']
   linoth[*,5] = ['[SII]6731','']
   argspltlin1 = {nx: 3, ny: 2,$
      label: ['[OII]3729','[NeIII]3869','Hbeta',$
      '[OI]6300','[NII]6548','[SII]6716'],$
      wave: [3726,3869,4950,6300,6563,6725],$
      off: [[-100,95],[-20,120],[-150,100],$
      [-100,75],[-80,100],[-100,80]],$
      linoth: linoth}
   ; Start with ap1 (not 3" aperture spectrum)
   for i=0,naps-1 do begin
      restore,'/Users/drupke/specfits/esi/makani/makani_00'+$
         string(i+1+apoff,format='(I02)')+'.xdr'
      emlz_ap = hash()
      foreach line, emlz->keys() do $
         emlz_ap[line] = reform(emlz[line,i+apoff,*,*],2)
      ;emlsig_ap = hash((emlsig->keys())[0],reform(emlsig['[OII]3729',i+apoff,*,*],2))
      ifsf_pltlin,struct,argspltlin1,plotdir+'makani-spec-'+$
         string(i+1,format='(I02)'),emlz_ap,/ps,$
         velinset=list('[OII]3729',[0.05d,0.53d,0.35d,0.83d],[-1499.99,999.99]),$
         yranscl=0.25d,boxsmooth=3d,title='Aperture '+string(i+1,format='(I0)')
   endfor

   ; --------------
   ; Aperture table for paper
   ; --------------

   openw,tab1,'/Users/drupke/ifs/esi/docs/makani_aptab.tex',/get_lun
   amp = ' & '
   dslash = ' \\'
   apsize = '1\arcsec$\times$' + $
      ['1\farcs54',$
      '1\farcs54',$
      '1\farcs54',$
      '1\farcs54',$
      '1\farcs54',$
      '3\farcs08',$
      '3\farcs08',$
      '1\farcs54',$
      '3\farcs08']
   for i=0,naps-1 do begin
      if apradkpc[i] eq 0 then apradkpc_lo_use = '' $
         else apradkpc_lo_use = string(apradkpc_lo[i],format='(D0.1)')
      printf,tab1,$
         'ap'+string(i+1,format='(I0)'),amp,apsize[i],amp,$
         apradkpc[i],'$_{',apradkpc_lo_use,'}^{',apradkpc_hi[i],'}$'+amp,$
         fx['[OII]3726+[OII]3729',i,2],'$\pm$',fxe['[OII]3726+[OII]3729',i,2],dslash,$
         format='(A-5,A3,A27,A3'+$
         ',D0.1,A3,A0,A3,D0.1,A5'+$
         ',D0.2,A5,D0.2,A3'+$
         ')'
   endfor
   
   free_lun,tab1
   
   ; --------------
   ; Flux table for paper
   ; --------------

   openw,tab3,'/Users/drupke/ifs/esi/docs/makani_fluxtab.tex',/get_lun
   amp = ' & '
   dslash = ' \\'
   tablines = ['[NeV]3426','[OII]3726','[OII]3729','[OII]3726+[OII]3729',$
      '[NeIII]3869','[OIII]4363','Hbeta','[OIII]5007',$
      '[NI]5198+[NI]5200','[NII]5755','[OI]6300','[NII]6583',$
      '[SII]6716','[SII]6731','[SII]6716+[SII]6731']
   ; if set to same as sigcut in IFSFA, just duplicates that and has no effect
   sigcut = 2d
   ; errors in log here are from d(log f1 - log f2) = sqrt((df1/f1)^2d+(df2/f2)^2)/ln 10
   ; d(log x) = 1/(x log 10), and d(g(f1,f2)) = sqrt((dg/df1 * df1)^2d + ...)
   ; observed fluxes
   for i=0,n_elements(tablines)-1 do begin
      drat1 = sqrt((fxe[tablines[i],0,0]/fx[tablines[i],0,0])^2d + $
         (fxe['Halpha',0,0]/fx['Halpha',0,0])^2d) / alog(10d)
      drat2 = sqrt((fxe[tablines[i],0,1]/fx[tablines[i],0,1])^2d + $
         (fxe['Halpha',0,1]/fx['Halpha',0,1])^2d) / alog(10d)
      drat = sqrt((fxe[tablines[i],0,2]/fx[tablines[i],0,2])^2d + $
         (fxe['Halpha',0,2]/fx['Halpha',0,2])^2d) / alog(10d)
      if fx[tablines[i],0,0]/fxe[tablines[i],0,0] gt sigcut then $
         c1str = string(alog10(fx[tablines[i],0,0]/fx['Halpha',0,0]),$
         '$\pm$',drat1,format='(D7.2,A5,D4.2)') $
      else $
         c1str = '\nodata'
      if fx[tablines[i],0,1]/fxe[tablines[i],0,1] gt sigcut then $
         c2str = string(alog10(fx[tablines[i],0,1]/fx['Halpha',0,1]),$
         '$\pm$',drat2,format='(D7.2,A5,D4.2)') $
      else $
         c2str = '\nodata'
      if fx[tablines[i],0,2]/fxe[tablines[i],0,2] gt sigcut then $
         ctotstr = string(alog10(fx[tablines[i],0,2]/fx['Halpha',0,2]),$
         '$\pm$',drat,format='(D7.2,A5,D4.2)') $
      else $
         ctotstr = '\nodata'
      tabstr = string(tablines[i],amp,c1str,amp,c2str,amp,ctotstr,$
         format='(A20,A3,A16,A3,A16,A3,A16)')
      for j=1,8 do begin
         if fx[tablines[i],j,2] ne bad AND $
            fx[tablines[i],j,2] ne 0d AND $
            fx['Halpha',j,2] ne bad AND $
            fx['Halpha',j,2] ne 0 then begin
            if fx[tablines[i],j,2]/fxe[tablines[i],j,2] gt sigcut then begin
               drat = sqrt((fxe[tablines[i],j,2]/fx[tablines[i],j,2])^2d + $
                  (fxe['Halpha',j,2]/fx['Halpha',j,2])^2d) / alog(10d)
               tabstr += string($
                  amp,alog10(fx[tablines[i],j,2]/fx['Halpha',j,2]),$
                  '$\pm$',drat,format='(A3,D7.2,A5,D4.2)')
            endif else begin
               tabstr += string(amp,'\nodata',format='(A3,A16)')
            endelse
         endif else begin
            tabstr += string(amp,'\nodata',format='(A3,A16)')
         endelse
      endfor
      if i ne n_elements(tablines)-1 then tabstr += dslash+'\relax' $
      else tabstr += dslash
      printf,tab3,tabstr
   endfor
   printf,tab3,'\hline'
   drat1 = fxe['Halpha',0,0]/fx['Halpha',0,0] / alog(10d)
   drat2 = fxe['Halpha',0,1]/fx['Halpha',0,1] / alog(10d)
   drat = fxe['Halpha',0,2]/fx['Halpha',0,2] / alog(10d)
   tabstr = string('log[f(Halpha)]',amp,$
      alog10(fx['Halpha',0,0])-16d,'$\pm$',drat1,amp,$
      alog10(fx['Halpha',0,1])-16d,'$\pm$',drat2,amp,$
      alog10(fx['Halpha',0,2])-16d,'$\pm$',drat,$
      format='(A30,A3,D7.3,A5,D5.3,A3,D7.3,A5,D5.3,A3,D7.3,A5,D5.3)')
   for j=1,8 do begin
      drat = fxe['Halpha',j,2]/fx['Halpha',j,2] / alog(10d)
      tabstr += string(amp,alog10(fx['Halpha',j,2])-16d,'$\pm$',drat,$
         format='(A3,D6.2,A5,D4.2)')
   endfor
   tabstr += dslash
   printf,tab3,tabstr
   printf,tab3,'\hline'

   ; intrinsic fluxes
   for i=0,n_elements(tablines)-1 do begin
      drat1 = sqrt((ecfxe[tablines[i],0,0]/ecfx[tablines[i],0,0])^2d + $
         (ecfxe['Halpha',0,0]/ecfx['Halpha',0,0])^2d) / alog(10d)
      drat2 = sqrt((ecfxe[tablines[i],0,1]/ecfx[tablines[i],0,1])^2d + $
         (ecfxe['Halpha',0,1]/ecfx['Halpha',0,1])^2d) / alog(10d)
      drat = sqrt((ecfxe[tablines[i],0,2]/ecfx[tablines[i],0,2])^2d + $
         (ecfxe['Halpha',0,2]/ecfx['Halpha',0,2])^2d) / alog(10d)
      if fx[tablines[i],0,0]/fxe[tablines[i],0,0] gt sigcut then $
         c1str = string(alog10(ecfx[tablines[i],0,0]/ecfx['Halpha',0,0]),$
         '$\pm$',drat1,format='(D7.2,A5,D4.2)') $
      else $
         c1str = '\nodata'
      if fx[tablines[i],0,1]/fxe[tablines[i],0,1] gt sigcut then $
         c2str = string(alog10(ecfx[tablines[i],0,1]/ecfx['Halpha',0,1]),$
         '$\pm$',drat2,format='(D7.2,A5,D4.2)') $
      else $
         c2str = '\nodata'
      if fx[tablines[i],0,2]/fxe[tablines[i],0,2] gt sigcut then $
         ctotstr = string(alog10(ecfx[tablines[i],0,2]/ecfx['Halpha',0,2]),$
         '$\pm$',drat,format='(D7.2,A5,D4.2)') $
      else $
         ctotstr = '\nodata'
      tabstr = string(tablines[i],amp,c1str,amp,c2str,amp,ctotstr,$
         format='(A20,A3,A16,A3,A16,A3,A16)')
      for j=1,8 do begin
         if ecfx[tablines[i],j,2] ne bad AND $
            ecfx[tablines[i],j,2] ne 0d AND $
            ecfx['Halpha',j,2] ne bad AND $
            ecfx['Halpha',j,2] ne 0 then begin
            if fx[tablines[i],j,2]/fxe[tablines[i],j,2] gt sigcut then begin
               drat = sqrt((ecfxe[tablines[i],j,2]/ecfx[tablines[i],j,2])^2d + $
                  (ecfxe['Halpha',j,2]/ecfx['Halpha',j,2])^2d) / alog(10d)
               tabstr += string($
                  amp,alog10(ecfx[tablines[i],j,2]/ecfx['Halpha',j,2]),$
                  '$\pm$',drat,format='(A3,D7.2,A5,D4.2)')
            endif else begin
               tabstr += string(amp,'\nodata',format='(A3,A16)')
            endelse
         endif else begin
            tabstr += string(amp,'\nodata',format='(A3,A16)')
         endelse
      endfor
      if i ne n_elements(tablines)-1 then tabstr += dslash+'\relax' $
      else tabstr += dslash
      printf,tab3,tabstr
   endfor
   printf,tab3,'\hline'
   tabstr = string('log[L(Halpha)]',amp,$
      loghalum[0,0],'$^{+',loghalumerrhi[0,0],'}_{-',loghalumerrlo[0,0],'}$',amp,$
      loghalum[0,1],'$^{+',loghalumerrhi[0,1],'}_{-',loghalumerrlo[0,1],'}$',amp,$
      loghalum[0,2],'$^{+',loghalumerrhi[0,2],'}_{-',loghalumerrlo[0,2],'}$',$
      format='(A20,A3,D5.2,A4,D4.2,A4,D4.2,A2,A3,D5.2,A4,D4.2,A4,D4.2,A2,A3,D5.2,A4,D4.2,A4,D4.2,A2)')
   for j=1,8 do begin
      tabstr += $
         string(amp,loghalum[j,2],'$^{+',loghalumerrhi[j,2],'}_{-',loghalumerrlo[j,2],'}$',$
         format='(A3,D5.2,A4,D4.2,A4,D4.2,A2)')
   endfor
   ;tabstr += dslash+'\relax'
   printf,tab3,tabstr


   free_lun,tab3

   ;----------------------------------------
   ;  Outflow rates from Voronoi-binned data
   ;----------------------------------------

   line = '[OII]3726+[OII]3729'
   linelab = ifsf_linesyntax(line)
   ; restore Voronoi fit; relaxed is for backwards compatibility with new version of IDL
   restore,file='/Users/drupke/ifs/kcwi/maps/hizeaj2118/vor/hizeaj2118.xdr',/relaxed
   ; start by separating into Episode I and II components; same def as in 2019 paper
   mapsig = windstr.e_vel['vsig',line]
   mapv50 = windstr.e_vel['v%50',line]
   mapv98 = windstr.e_vel['v%98',line]
   igd = where(mapv98 ne bad,ctgd)
   iep1 = where(mapv98 gt -700 AND mapv98 lt bad)
   xind = rebin(indgen(kcwidx)+1,kcwidx,kcwidy)
   iep2 = where(mapv98 le -700 AND xind gt 15)
   ; Observed [OII] fluxes in 10^-16 erg/s/cm^-2, per pixel
   mapvorf = windstr.e_flx['ftot',line]
   igdvorf = where(mapvorf ne bad)
   mapvorf[igdvorf] = mapvorf[igdvorf]*1d-16
   ; Placeholder for Halpha intrinsic fluxes
   mapvorfha = mapvorf[*]
   mapvorfha_lo = mapvorf[*]
   mapvorfha_hi = mapvorf[*]
   ; x, y arrays in 0-offset coordinates
   map_x = rebin(dindgen(kcwidx)+1,kcwidx,kcwidy)
   map_y = rebin(transpose(dindgen(kcwidy)+1),kcwidx,kcwidy)
   ; pixel radius arrays 
   map_r = sqrt((map_x - kcwi_center_nuclei[0])^2d + (map_y - kcwi_center_nuclei[1])^2d)


   ; Voronoi coordinates
   initdatvor=call_function('ifsf_hizeaj2118vor')
   vormap=initdatvor.vormap
   nvorind = max(vormap)
   vorcoords = intarr(nvorind,2)
   ; average radius for each Voronoi bin in kpc
   rvor = dblarr(kcwidx,kcwidy)+bad
   for i=1,nvorind do begin
      ivor = where(vormap eq i,ctivor)
      rvor[ivor] = mean(map_r[ivor])*kcwi_kpc_per_pix
      ; Use [OII]/Halpha fit to compute Halpha
      ;if rvor[ivor[0]] le xfitlo then $
      ;   o2haconv = 10d^yfitlo $
      ;else if rvor[ivor[0]] gt xfitlo AND rvor[ivor[0]] lt xfithi then $
      ;   o2haconv = 10d^(o2ha_exto2_pars[0]+rvor[ivor[0]]*o2ha_exto2_pars[1]) $
      ;else $
      ;    o2haconv = 10d^yfithi
      o2haconv = 10d^(o2ha_exto2_pars[0]+rvor[ivor[0]]^fitpower*o2ha_exto2_pars[1])
      ;print,i,rvor[ivor[0]],o2haconv
      mapvorfha[ivor] /= o2haconv
      if rvor[ivor[0]] lt 5d then begin
         ; correct fluxes inside 5kpc for contribution from star formation
         mapvorfha[ivor]*=0.6d
         ; inside this radius, assume model is perfectly correct
         mapvorfha_lo[ivor] /= o2haconv
         mapvorfha_hi[ivor] /= o2haconv
      endif else begin
         ; outside 5kpc, estimate model error by adding in RMS
         mapvorfha_hi[ivor] /= o2haconv/2d
         mapvorfha_lo[ivor] /= o2haconv*2d
      endelse
   endfor

   ; Use [OII]/Halpha fit to compute Halpha

   lum = hash('ep1',0d,'ep2',0d)
   mass = hash('ep1',0d,'ep2',0d)
   vrad =  hash('ep1',0d,'ep2',0d)
   massdot = hash('ep1',0d,'ep2',0d)
   mom = hash('ep1',0d,'ep2',0d)
   momdot = hash('ep1',0d,'ep2',0d)
   ener = hash('ep1',0d,'ep2',0d)
   enerdot = hash('ep1',0d,'ep2',0d)
   vrad98 =  hash('ep1',0d,'ep2',0d)
   massdot98 = hash('ep1',0d,'ep2',0d)
   mom98 = hash('ep1',0d,'ep2',0d)
   momdot98 = hash('ep1',0d,'ep2',0d)
   ener98 = hash('ep1',0d,'ep2',0d)
   enerdot98 = hash('ep1',0d,'ep2',0d)

   map_lum = hash('ep1',dblarr(kcwidx,kcwidy),'ep2',dblarr(kcwidx,kcwidy))
   map_mass = hash('ep1',dblarr(kcwidx,kcwidy),'ep2',dblarr(kcwidx,kcwidy))
   map_vrad = hash('ep1',dblarr(kcwidx,kcwidy),'ep2',dblarr(kcwidx,kcwidy))
   map_mom = hash('ep1',dblarr(kcwidx,kcwidy),'ep2',dblarr(kcwidx,kcwidy))
   map_ener = hash('ep1',dblarr(kcwidx,kcwidy),'ep2',dblarr(kcwidx,kcwidy))
   map_massdot = hash('ep1',dblarr(kcwidx,kcwidy),'ep2',dblarr(kcwidx,kcwidy))
   map_momdot = hash('ep1',dblarr(kcwidx,kcwidy),'ep2',dblarr(kcwidx,kcwidy))
   map_enerdot = hash('ep1',dblarr(kcwidx,kcwidy),'ep2',dblarr(kcwidx,kcwidy))
   map_vrad98 = hash('ep1',dblarr(kcwidx,kcwidy),'ep2',dblarr(kcwidx,kcwidy))
   map_mom98 = hash('ep1',dblarr(kcwidx,kcwidy),'ep2',dblarr(kcwidx,kcwidy))
   map_ener98 = hash('ep1',dblarr(kcwidx,kcwidy),'ep2',dblarr(kcwidx,kcwidy))
   map_massdot98 = hash('ep1',dblarr(kcwidx,kcwidy),'ep2',dblarr(kcwidx,kcwidy))
   map_momdot98 = hash('ep1',dblarr(kcwidx,kcwidy),'ep2',dblarr(kcwidx,kcwidy))
   map_enerdot98 = hash('ep1',dblarr(kcwidx,kcwidy),'ep2',dblarr(kcwidx,kcwidy))
   iepgd = hash('ep1',iep1,'ep2',iep2)
   
   ; wind radii in kpc
   starage = hash('ep1',4d8,'ep2',7d6)
   radkpc = hash('ep1',40d,'ep2',15d)
   ;elecden = hash('ep1',dblarr(kcwidx,kcwidy)+1d,'ep2',dblarr(kcwidx,kcwidy)+1000d)
   elecden_pre = hash('ep1',dblarr(kcwidx,kcwidy)+1d,'ep2',dblarr(kcwidx,kcwidy)+10d)
   elecden_post = hash('ep1',dblarr(kcwidx,kcwidy)+100d,'ep2',dblarr(kcwidx,kcwidy)+1000d)

   ; velocities in cm/s
   mapsig *= 1d5
   mapv50 *= 1d5
   mapv98 *= 1d5

   foreach key, radkpc->keys() do begin
      ; ... in cm
      radcm = 3.08567802d18 * radkpc[key] * 1d3
      ; ... in spaxels
      radpix = (radkpc[key] / kpc_per_as)/kcwips
      iiepgd = cgsetintersection(where(map_r+0.5d lt radpix),iepgd[key])
      ; Extra +0.5 is to prevent map_dtheta from being undefined, and will
      ; introduce an error of only half a pixel. OK if Rpix >> 1.
      map_dtheta = asin((map_r+0.5d)/radpix) - asin((map_r-0.5d)/radpix)
      map_dphi_sintheta = 1d/radpix
      map_domega = map_dphi_sintheta*map_dtheta
      map_costheta = cos(asin(map_r/radpix))

      map_vrad[key,iiepgd] = abs(mapv50[iiepgd]) / map_costheta[iiepgd]
      map_vrad98[key,iiepgd] = abs(mapv98[iiepgd]) / map_costheta[iiepgd]

      map_lum[key,iiepgd] = drt_linelum(mapvorfha[iiepgd]*1d-3,dist,/ergs)
      map_mass[key,iiepgd] = $
         drt_hiimass(map_lum[key,iiepgd]/2d,elecden_pre[key,iiepgd]) + $
         drt_hiimass(map_lum[key,iiepgd]/2d,elecden_post[key,iiepgd])       
      ;map_mass[key,iiepgd] = drt_hiimass(map_lum[key,iiepgd],elecden[key,iiepgd])
      ; The following formula comes from dividing dM/dt^avg_thin by M_thin
      ; (eq. 7 and 8 in RVS05b). Note that this is basically equivalent to
      ; dividing each pixel by its own dynamical time, tdyn ~ R/v.
      map_massdot[key,iiepgd] = map_mass[key,iiepgd] * map_vrad[key,iiepgd] * $
         speryr / radcm
      ; in units of msun/yr
      map_mom[key,iiepgd] = map_mass[key,iiepgd] * msun * map_vrad[key,iiepgd]
      ; in units of dyne s
      map_momdot[key,iiepgd] = map_massdot[key,iiepgd] * msun / speryr * $
         map_vrad[key,iiepgd]
      ; dp/dt, in Lsun
      map_ener[key,iiepgd] = map_mass[key,iiepgd] * msun * $
         (0.5d*map_vrad[key,iiepgd]^2d + 1.5d*mapsig[iiepgd]^2d)
      ; in erg
      map_enerdot[key,iiepgd] = map_massdot[key,iiepgd] * msun / speryr * $
         (0.5d*map_vrad[key,iiepgd]^2d + 1.5d*mapsig[iiepgd]^2d)
      ; in erg/s
      map_massdot98[key,iiepgd] = map_mass[key,iiepgd] * map_vrad98[key,iiepgd] * $
         speryr / radcm
      ; in units of msun/yr
      map_mom98[key,iiepgd] = map_mass[key,iiepgd] * msun * map_vrad98[key,iiepgd]
      ; in units of dyne s
      map_momdot98[key,iiepgd] = map_massdot98[key,iiepgd] * msun / speryr * $
         map_vrad98[key,iiepgd]
      ; dp/dt*c, in Lsun
      map_ener98[key,iiepgd] = map_mass[key,iiepgd] * msun * $
         (0.5d*map_vrad98[key,iiepgd]^2d + 1.5d*mapsig[iiepgd]^2d)
      ; in erg
      map_enerdot98[key,iiepgd] = map_massdot98[key,iiepgd] * msun / speryr * $
         (0.5d*map_vrad98[key,iiepgd]^2d + 1.5d*mapsig[iiepgd]^2d)
      ; in erg/s
      
      lum[key] = total(map_lum[key,iiepgd])
      mass[key] = total(map_mass[key,iiepgd])
      massdot[key] = total(map_massdot[key,iiepgd])
      mom[key] = total(map_mom[key,iiepgd])
      momdot[key] = total(map_momdot[key,iiepgd])
      ener[key] = total(map_ener[key,iiepgd])
      enerdot[key] = total(map_enerdot[key,iiepgd])
      massdot98[key] = total(map_massdot98[key,iiepgd])
      mom98[key] = total(map_mom98[key,iiepgd])
      momdot98[key] = total(map_momdot98[key,iiepgd])
      ener98[key] = total(map_ener98[key,iiepgd])
      enerdot98[key] = total(map_enerdot98[key,iiepgd])

      ;print,key,lum[key],mass[key],massdot[key],mom[key],momdot[key],ener[key],enerdot[key]
;      print,'mean deprojected velocity from v50 and v98'
;      print,key,mean(map_vrad[key,iiepgd])/1d5
;      print,key,mean(map_vrad98[key,iiepgd])/1d5
      vrad[key] = mean(map_vrad[key,iiepgd])/1d5
      vrad98[key] = mean(map_vrad98[key,iiepgd])/1d5
;
;      print,'mdot from dynamical time, corresonding flow velocity'
;      print,key,mass[key]/starage[key],$
;         radkpc[key]*1d3*3.08d13/(starage[key]*365.25*24d*3600d)

   endforeach

   ; total Ha luminosity
   lumtot_wind = total(drt_linelum(mapvorfha[igdvorf]*1d-3,dist,/ergs))
   lumtot_wind_lo = total(drt_linelum(mapvorfha_lo[igdvorf]*1d-3,dist,/ergs))
   lumtot_wind_hi = total(drt_linelum(mapvorfha_hi[igdvorf]*1d-3,dist,/ergs))
   print,'Total Halpha luminosity:',lumtot_wind,' + ',lumtot_wind_hi-lumtot_wind,' - ',lumtot_wind-lumtot_wind_lo

   openw,tab2,'/Users/drupke/ifs/esi/docs/makani_masstab.tex',/get_lun
   amp = ' & '
   dslash = ' \\'
   nodat = ' \nodata'
   epsrom = ['I','II']
   for i=0,1 do begin
      ep = 'ep'+string(i+1,format='(I0)')
      meth = '1a'
      ;if i eq 1 then meth = '1a'
      printf,tab2,epsrom[i],amp,meth,amp,lum[ep],amp,mass[ep],amp,$
         vrad[ep],amp,massdot[ep],$
         amp,mom[ep],amp,momdot[ep],amp,ener[ep],amp,enerdot[ep],dslash,$
         format='(A-3,A3,A3,A3,E9.2,A3,E9.2,A3,I5,A3,D6.1,A3,'+$
         'E9.2,A3,E9.2,A3,E9.2,A3,E9.2,A3)'
      ;if i eq 1 then begin
         meth = '1b'
         printf,tab2,epsrom[i],amp,meth,amp,lum[ep],amp,mass[ep],amp,$
            vrad98[ep],amp,massdot98[ep],$
            amp,mom98[ep],amp,momdot98[ep],amp,ener98[ep],amp,enerdot98[ep],dslash,$
            format='(A-3,A3,A3,A3,E9.2,A3,E9.2,A3,I5,A3,D6.1,A3,'+$
            'E9.2,A3,E9.2,A3,E9.2,A3,E9.2,A3)'
      ;endif
      meth = '2'
      printf,tab2,epsrom[i],amp,meth,amp,lum[ep],amp,mass[ep],amp,$
         radkpc[ep]*1d3*3.08d13/(starage[ep]*365.25*24d*3600d),amp,$
         mass[ep]/starage[ep],$
         amp,nodat,amp,nodat,amp,nodat,amp,nodat,dslash,$
         format='(A-3,A3,A3,A3,E9.2,A3,E9.2,A3,I5,A3,D6.1,A3,'+$
         'A9,A3,A9,A3,A9,A3,A9,A3)'
      if i eq 0 then printf,tab2,'\hline'
   endfor
   
   free_lun,tab2
   
end