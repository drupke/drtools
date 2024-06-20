; docformat = 'rst'
;
;+
;
; Process data for ORCs spectra paper.
;
; :Categories:
;    DRTOOLS
;
; :Returns:
;    None.
;
; :Params:
;
; :Keywords:
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
;      2023, DSNR, created
;
; :Copyright:
;    Copyright (C) 2023-2024 David S. N. Rupke
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
pro orcs_plots

   bad = 1d99
   
   speryr = 24d*3600d*365.25d    ; seconds in a year
   lsun = 3.826d33               ; solar luminosities, erg/s
   msun = 1.989e33               ; solar mass, g
   c_cms = 2.99792d10

   plotdir = '/Users/drupke/Box Sync/orc-spectra/plots/'
   plotquantum = 2.5 ; in inches

   ; Planck cosmology
   angdist=1b
   dist = drt_plumdist(0.451d,angdist=angdist)

   kpc_per_as = angdist*1000d/206265d
   gmosps = 0.1614d
   kpc_per_pix = gmosps * kpc_per_as

   ; reddening
   rv = 4.05d

   ; With 4se
;   orclabs = ['1c','23c','4c','4se','5c','5se']
;   spectype = ['LINER','AGN','LINER','\nodata','HII/LINER','HII']
;   orccol = ['black','cyan','blue','red','magenta','brown']
;   tabfile = '/Users/drupke/Box Sync/orc-spectra/tabs/tab-spec.tex'
;   plotsuf = '-with4se'
   ; ... and without
   orclabs = ['1c','23c','4c','5c','5se']
   spectype = ['LINER','AGN','LINER','HII/LINER','HII']
   orccol = ['black','cyan','blue','red','magenta']
   orcsymsize = [1.75d,1d,1.75d,1.75d,1d]
   tabfile = '/Users/drupke/Box Sync/orc-spectra/tabs/tab-spec-no4se.tex'
   plotsuf = ''
   ; apertures
   naps = n_elements(orclabs)
   ; Get data
   linelab = 1b
   linelist = ifsf_linelist(/all,linelab=linelab,/vacuum)

   ; fluxes and errors of lines
   fx = hash()
   fxe = hash()
   weq = hash()
   weqe = hash()
   ecfx = hash()
   ecfxe = hash()
   ; gas zs and sigmas
   z = dblarr(naps) + bad
   z_err = dblarr(naps) + bad
   sig = dblarr(naps) + bad
   sig_err = dblarr(naps) + bad
   ; populate stellar properties
   stel_z = dblarr(naps) + bad
   stel_z_err = dblarr(naps) + bad
   stel_sigma = dblarr(naps) + bad
   stel_sigma_err = dblarr(naps) + bad
   stel_ebv = dblarr(naps) + bad
   stel_ebv_err = dblarr(naps) + bad
   stel_percents = dblarr(naps,3) + bad
   ; line ratios
   lrs = hash()
   lrsel = hash()
   lrseh = hash()
   lrselin = hash()
   eclrs = hash()
   eclrsel = hash()
   eclrseh = hash()
   eclrselin = hash()
   ;
   dn4000 = dblarr(naps) + 0d
   dn4000err = dblarr(naps) + 0d
   ; densities
   dens_s2 = dblarr(naps) + bad
   dens_s2errlo = dblarr(naps) + bad
   dens_s2errhi = dblarr(naps) + bad
   dens_o2 = dblarr(naps) + bad
   dens_o2errlo = dblarr(naps) + bad
   dens_o2errhi = dblarr(naps) + bad
   ; Halpha luminosities
   loghalum = dblarr(naps) + bad
   loghalumerrlo = dblarr(naps) + bad
   loghalumerrhi = dblarr(naps) + bad
 
   ncomp=1
   for iap=0,naps-1 do begin
      ; fluxes, line ratios
      restore,'/Users/drupke/specfits/gmos/orc'+orclabs[iap]+$
         '/orc'+orclabs[iap]+'.lin.xdr'
      ; redshifts
      restore,'/Users/drupke/specfits/gmos/orc'+orclabs[iap]+$
         '/orc'+orclabs[iap]+'.lininit.xdr'
      ; stellar properties
      restore,'/Users/drupke/specfits/gmos/orc'+orclabs[iap]+$
         '/orc'+orclabs[iap]+'.cont.xdr'
      ; Compute CVDF velocities
      emlcvdfvel = ifsf_cmpcvdfvals(emlcvdf,emlflx,emlflxerr)
      ; populate fluxes and errors of lines
      if iap eq 0 then begin
         foreach line, linelist.keys() do begin
            fx[line] = dblarr(naps)+bad
            fxe[line] = dblarr(naps)+bad
            weq[line] = dblarr(naps)+bad
         endforeach
      endif
      foreach key,emlflx['fc1'].keys() do begin
         fx[key,iap] = emlflx['fc1',key,0]
         fxe[key,iap] = emlflxerr['fc1',key,0]
         weq[key,iap] = emlweq['fc1',key,0]
      endforeach
      ; For 5c, set undetected [OI] equal to [SII]6731 detection as a conservative
      ; upper limit. Set error to value.
      if orclabs[iap] eq '5c' then begin
         fx['[OI]6300',iap] = emlflx['fc1','[SII]6731',0]
         fxe['[OI]6300',iap] = emlflx['fc1','[SII]6731',0]
      endif
      ; populate gas zs and sigmas
      ;z_tmp = reform(emlz['[OII]3729',iapp1,*,*])
      if emlwav['c1','Halpha',0] ne bad then begin
         z[iap] = (emlwav['c1','Halpha',0]/linelist['Halpha']) - 1d
         z_err[iap] = (emlwaverr['c1','Halpha',0]/linelist['Halpha'])
         sig[iap] = emlsig['c1','Halpha',0]
         sig_err[iap] = emlsigerr['c1','Halpha',0]
      endif
;     populate stellar properties
      stel_z[iap] = contcube.stel_z[0]
      stel_z_err[iap] = mean(contcube.stel_z_err)
      stel_sigma[iap] = contcube.stel_sigma[0]
      stel_sigma_err[iap] = mean(contcube.stel_sigma_err)
      stel_ebv[iap] = contcube.stel_ebv[0]
      ; put lower limit of 0.01 on E(B-V) error
      stel_ebv_err[iap] = (mean(contcube.stel_ebv_err) le 0.01) OR $
          (mean(contcube.stel_ebv_err) eq bad) ? 0.01 : mean(contcube.stel_ebv_err)
      stel_percents[iap,*] = $
         [contcube.stel_percents[2]+contcube.stel_percents[3],$
          contcube.stel_percents[4],contcube.stel_percents[5]]

      ; Compute dn4000
      restore,'/Users/drupke/specfits/gmos/orc'+orclabs[iap]+$
         '/orc'+orclabs[iap]+'_0001.xdr'
      zwave = struct.wave/(1d +stel_z[iap])
      ibluwin = where(zwave ge 3850d AND zwave le 3950d,ctblu)
      iredwin = where(zwave ge 4000d AND zwave le 4100d,ctred)
      if ctred gt 0 AND ctblu gt 0 then begin
         redmean = mean(struct.spec[iredwin])
         ;redmeanerr = sqrt(total(struct.spec_err[iredwin]^2d))/ctred
         redmeanerr = stddev(struct.spec[iredwin])/sqrt(ctred)
         blumean = mean(struct.spec[ibluwin])
         ;blumeanerr = sqrt(total(struct.spec_err[ibluwin]^2d))/ctblu
         blumeanerr = stddev(struct.spec[ibluwin])/sqrt(ctblu)
         dn4000[iap] = redmean/blumean
         dn4000err[iap] = sqrt(redmeanerr^2d + blumeanerr^2d)
       endif
   endfor

   ; compute extinction
   ; First do sigma cut on all lines
   fx_use = fx[*]
   fxe_use = fxe[*]
   foreach key,fx.keys() do begin
      ilowsig = where(fx_use[key] le 2.5*fxe_use[key] AND fx_use[key] ne bad,ctlowsig)
      if ctlowsig gt 0 then begin
         fx_use[key,ilowsig] = bad
         fxe_use[key,ilowsig] = bad
      endif
   endforeach
   lrs_ebv = ifsf_lineratios(fx_use,fxe_use,linelist,errlo=lrsel,errhi=lrseh,$
      errlin=lrselin,/ebvonly,rv=rv)
  
   ; Now extincted line ratio, without 1sigma cut
   lrs = ifsf_lineratios(fx,fxe,linelist,errlo=lrsel,errhi=lrseh,$
      errlin=lrselin,rv=rv)
  
   for iap=0,naps-1 do begin
      ebv_use = lrs_ebv['ebv',iap]
      ebverr_use = lrselin['ebv',iap]
      ; Substitute stellar E(B-V) if no gas measurement
      if lrs_ebv['ebv',iap] eq bad then begin
         ebv_use = stel_ebv[iap]/0.44
         ebverr_use = stel_ebv_err[iap]/0.44
      endif
      ; apply E(B-V)
      foreach key,fx.keys() do begin
         if iap eq 0 then begin
            ecfx[key] = dblarr(naps)+bad
            ecfxe[key] = dblarr(naps)+bad
         endif
         igd = where(fx[key,iap] ne bad AND ebv_use ne bad,ctgd)
         if ctgd gt 0 then begin
            ec_tmp = drt_dustcor_ccm(linelist[key],fx[key,iap],$
               ebv_use,ebverr=ebverr_use,$
               fluxerr=fxe[key,iap],rv=rv)
            ecfx[key,iap] = ec_tmp[*,0]
            ecfxe[key,iap] = ec_tmp[*,1]
         endif
      endforeach
      ; Use Halpha instead of Hbeta for 5c, since it is more significant.
      ; Stellar extinction is 0.
      ; Result is basically the same; point goes down with Halpha by about 0.2-0.3dex
      if orclabs[iap] eq '5c' then begin
         ecfx['Hbeta',iap] = ecfx['Halpha',iap] / 2.86
         ecfxe['Hbeta',iap] = ecfxe['Halpha',iap] / 2.86
      endif
   endfor
   
   ; compute unextincted line ratios
   lrlabs = ['n2ha','o1ha','s2ha','o3hb','o3o2','s2','o2','n1hb']
   lrlines = [['[NII]6583','Halpha'],$
      ['[OI]6300','Halpha'],$
      ['[SII]6716+[SII]6731','Halpha'],$
      ['[OIII]5007','Hbeta'],$
      ['[OIII]5007','[OII]3726+[OII]3729'],$
      ['[SII]6716','[SII]6731'],$
      ['[OII]3729','[OII]3726'],$
      ['[NI]5198+[NI]5200','Hbeta']]
   lrlist = hash()
   for i=0,n_elements(lrlabs)-1 do lrlist[lrlabs[i]]=lrlines[*,i]
   eclrs = ifsf_lineratios(ecfx,ecfxe,linelist,errlo=eclrsel,errhi=eclrseh,$
      errlin=eclrselin,/lronly,lrlist=lrlist,rv=rv)

   drt_voplot,lrs,plotdir+'orcs_vo_noextcor'+plotsuf+'.eps',errlo=lrsel,errhi=lrseh,$
      color=orccol,labels=orclabs,symsize=orcsymsize
   drt_voplot,eclrs,plotdir+'orcs_vo'+plotsuf+'.eps',errlo=eclrsel,errhi=eclrseh,$
      color=orccol,labels=orclabs,symsize=orcsymsize

   ; File for data-behind-the-figure for Figure 4
   fig4tab = '/Users/drupke/Box Sync/orc-spectra/tabs/tab-fig4.dat'
   openw,tablun,fig4tab,/get_lun
   printf,tablun,'# label','o3hb','siglo','sighi','n2ha','siglo','sighi',$
      's2ha','siglo','sighi','o1ha','siglo','sighi','o3o2','siglo','sighi',$
      format='(A-8,15A8)'
   tabrats = ['o3hb','n2ha','s2ha','o1ha','o3o2']
   tabbad = -9.9999
   for iap=0,naps-1 do begin
      eclrsel_o3hb = eclrsel['o3hb',iap]
      tabstr = string('ORC'+orclabs[iap],format='(A-8)')
      foreach trat,tabrats do begin
         thisrat = eclrs[trat,iap]
         thisratel = eclrsel[trat,iap]
         thisrateh = eclrseh[trat,iap]
         if thisrat eq bad then thisrat = tabbad
         if thisratel eq bad then thisratel = tabbad
         if thisrateh eq bad then thisrateh = tabbad
         tabstr += string(thisrat,thisratel,thisrateh,format='(3D8.4)')
      endforeach
      printf,tablun,tabstr
   endfor
   free_lun,tablun

;  Halpha luminosities
   hasfr = dblarr(naps) + bad
   hasfrerr = dblarr(naps) + bad
   loghalum = dblarr(naps) + bad
   loghalumerrlo = dblarr(naps) + bad
   loghalumerrhi = dblarr(naps) + bad
   igdha = where(ecfx['Halpha'] ne bad AND ecfxe['Halpha'] ne bad)
   haflux = ecfx['Halpha',igdha]*1d-17*1d-3 ; fluxes in W m^-2
   hafluxerr = ecfxe['Halpha',igdha]*1d-17*1d-3 ; fluxes in W m^-2
   halum = drt_linelum(haflux,dist,/ergs,err=hafluxerr,$
      lumerr=halumerr)
   ; https://ned.ipac.caltech.edu/level5/Sept12/Calzetti/Calzetti1_2.html
   hasfr[igdha] = 7.9d-42 * halum
   hasfrerr[igdha] = 7.9d-42 * halumerr
   loghalum[igdha] = alog10(halum)
   loghalumerrlo[igdha] = loghalum[igdha] - alog10(halum - halumerr)
   loghalumerrhi[igdha] = alog10(halum+halumerr)-loghalum[igdha]
   print,string(halum,format='(E0.2)')+'+/-'+$
      string(halumerr,format='(E0.2)')
;  o3 luminosities
   logo3lum = dblarr(naps) + bad
   logo3lumerrlo = dblarr(naps) + bad
   logo3lumerrhi = dblarr(naps) + bad
   igdo3 = where(ecfx['[OIII]5007'] ne bad AND ecfxe['[OIII]5007'] ne bad)
   o3flux = ecfx['[OIII]5007',igdo3]*1d-17*1d-3 ; fluxes in W m^-2
   o3fluxerr = ecfxe['[OIII]5007',igdo3]*1d-17*1d-3 ; fluxes in W m^-2
   o3lum = drt_linelum(o3flux,dist,/ergs,err=o3fluxerr,$
      lumerr=o3lumerr)
   logo3lum[igdo3] = alog10(o3lum)
   logo3lumerrlo[igdo3] = logo3lum[igdo3] - alog10(o3lum - o3lumerr)
   logo3lumerrhi[igdo3] = alog10(o3lum+o3lumerr)-logo3lum[igdo3]

   ; Black hole mass
   mbh = 8.32 + 5.64*alog10(stel_sigma/200.)
   mbherrlo = mbh - (8.32 + 5.64*alog10((stel_sigma - stel_sigma_err)/200.))
   mbherrup = (8.32 + 5.64*alog10((stel_sigma + stel_sigma_err)/200.)) - mbh
   mbherr = sqrt(0.3^2 + (mean([[mbherrlo],[mbherrup]],dim=2))^2)


; -----------
; Spec props table
; -----------

   openw,tabs,tabfile,/get_lun
   amp = ' & '
   dslash = ' \\'
   plmi = '\pm'
   stab = ''
   for iap=0, naps-1 do begin
      stab = ''
      stab += string('ORC'+orclabs[iap],format='(A8)')
      stab += string(amp,stel_z[iap],plmi,stel_z_err[iap],$
         format='(A3,D7.5,A3,D-7.5)')
      if z[iap] ne bad then $
        stab += string(amp,z[iap],plmi,z_err[iap],$
           format='(A3,D7.5,A3,D-7.5)') $
      else $
        stab += string(amp,'\nodata',format='(A3,A17)')
      stab += string(amp,stel_sigma[iap],plmi,stel_sigma_err[iap],$
         format='(A3,I3,A3,I-2)')
      if sig[iap] ne bad AND sig[iap] gt 50d then $
         stab += string(amp,sig[iap],plmi,sig_err[iap],$
            format='(A3,I3,A3,I-2)') $
      else $
         stab += string(amp,'\nodata',format='(A3,A8)')
      stab += string(amp,spectype[iap],format='(A3,A10)')
      if logo3lum[iap] ne bad then begin
         if finite(logo3lumerrlo[iap]) then $
            stab += $
            string(amp,logo3lum[iap],'^{+',logo3lumerrhi[iap],$
            '}_{-',logo3lumerrlo[iap],'}',$
            format='(A3,D5.2,A3,D4.2,A4,D4.2,A1)') $
         else $
            stab += $
            string(amp,'<',logo3lum[iap],'^{+',logo3lumerrhi[iap],'}',$
            format='(A3,A8,D-5.2,A3,D4.2,A-1)')
      endif else $
         stab += string(amp,'\nodata',format='(A3,A21)')
      if loghalum[iap] ne bad then $
         stab += $
            string(amp,loghalum[iap],'^{+',loghalumerrhi[iap],$
            '}_{-',loghalumerrlo[iap],'}',$
            format='(A3,D5.2,A3,D4.2,A4,D4.2,A1)') $
      else $
         stab += string(amp,'\nodata',format='(A3,A21)')
      if hasfr[iap] ne bad then begin
         if spectype[iap] eq 'HII' then $
            stab += $
            string(amp,hasfr[iap],plmi,hasfrerr[iap],$
            format='(A3,D5.1,A3,D-3.1)') $
         else $
            stab += $
            string(amp,'<',hasfr[iap],'^{+',loghalumerrhi[iap],'}',$
            format='(A3,A1,D-3.1,A3,D3.1,A1)')
      endif else $
         stab += string(amp,'\nodata',format='(A3,A11)')
      stab += string(amp,stel_ebv[iap],plmi,stel_ebv_err[iap],$
         format='(A3,D5.2,A3,D-4.2)')
      if lrs_ebv['ebv',iap] ne bad then $
         stab += string(amp,lrs_ebv['ebv',iap],plmi,lrselin['ebv',iap],$
            format='(A3,D5.2,A3,D-4.2)') $
      else $
         stab += string(amp,'\nodata',format='(A3,A12)')
      stab += string(amp,dn4000[iap],plmi,dn4000err[iap],format='(A3,D4.2,A3,D4.2)')
      stab += string(amp,stel_percents[iap,0],format='(A3,I3)')
      stab += string(amp,stel_percents[iap,1],format='(A3,I3)')
      stab += string(amp,stel_percents[iap,2],format='(A3,I3)')
      stab += string(amp,mbh[iap],plmi,mbherr[iap],$
         format='(A3,D5.2,A3,D-5.2)')
      stab += dslash
      printf,tabs,stab
   endfor
   free_lun,tabs
  
   for iap=0, naps-1 do begin
      print,orclabs[iap],weq['[OII]3726+[OII]3729',iap]/(1d + stel_z[iap]),$
         weq['[OIII]5007',iap]/(1d + stel_z[iap])
   endfor
   
end