; docformat = 'rst'
;
;+
;
; Plot MgII linemaps for a couple of Hi-z E+A galaxies.
;
; :Categories:
;    DRTOOLS/KCWI
;
; :Returns:
;    Postscript plots.
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
;      2019jun26, DSNR, created
;
; :Copyright:
;    Copyright (C) 2019 David S. N. Rupke
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
;
pro drt_hizeakcwimaps

  platescale = 0.2914d
  c_kms = 299792.458d
  linelistmg2 = ifsf_linelist(['MgII2795','MgII2802'])
  dir = '/Users/drupke/ifs/kcwi/'
  samplefac = 10d
  logconstant = 100d

; J0905

  zsys = 0.7114
  label = 'hizeaj0905'
  infile = dir+'cubes/'+label+'/'+label+'.fits'

  ; Planck 2018 parameters: https://ui.adsabs.harvard.edu/#abs/arXiv:1807.06209
  ldist = lumdist(zsys,H0=67.4d,Omega_m=0.315d,Lambda0=0.685d,/silent)
  asdist = ldist/(1d + zsys)^2d
  kpc_per_as = asdist*1000d/206265d
  kpc_per_pix = platescale * kpc_per_as


  cube = ifsf_readcube(infile,/quiet,header=header,datext=-1,varext=1,dqext=2)

  size_tmp = size(cube.dat)
  dx = size_tmp[1]
  dy = size_tmp[2]
  dz = size_tmp[3]
  center_axes = [double(dx)/2d,double(dy)/2d]+0.5d
  center_nuclei = [30d,30d]

;  coordinates in kpc
   xran_kpc = double([-(center_nuclei[0]-0.5),dx-(center_nuclei[0]-0.5)]) $
              * kpc_per_pix
   yran_kpc = double([-(center_nuclei[1]-0.5),dy-(center_nuclei[1]-0.5)]) $
              * kpc_per_pix
   center_nuclei_kpc_x = 0d
   center_nuclei_kpc_y = 0d

; For IDL-Colorbars to work properly
  device, true=24,decompose=0,retain=2
  ;  For interpolation + subsampling
  xresamp = (dindgen(dx*samplefac)-samplefac/2d + 0.5d)/samplefac
  yresamp = (dindgen(dy*samplefac)-samplefac/2d + 0.5d)/samplefac

  vlim = [-800,1200]
  redwave = (1d + zsys)*((linelistmg2['MgII2795']+linelistmg2['MgII2802'])/2d)
  wavesum = redwave + vlim/c_kms * redwave
  wavesub = [[wavesum[0]-10d,wavesum[0]-0d],$
             [wavesum[1]+0d,wavesum[1]+10d]]
  mg2map = ifsr_makelinemap(cube,wavesum,wavesub=wavesub,/allow)
  maxmg2map = max(mg2map.dat)
  nmg2map = mg2map.dat / maxmg2map
  mg2map_fine = interpolate(nmg2map,xresamp,yresamp,/double,/grid,cubic=-0.5)
  maxmg2map /= platescale^2d
  maxmg2map *= 1d-16

  cmap = ifsr_makelinemap(cube,[3500,5500])
  maxcmap = max(cmap.dat)
  ncmap = cmap.dat / maxcmap
  cmap_fine = interpolate(ncmap,xresamp,yresamp,/double,/grid,cubic=-0.5)
  maxcmap /= (0.2914d)^2d
  maxcmap *= 1d-16


  cgps_open,dir+'maps/'+label+'/rb1/'+label+'mg2.eps',charsize=1d,$
            /encap,/inches,xs=7.5d,ys=7.5d,/qui,/nomatch

  zran = [0.1,1]
  ncbdiv = 5
  dzran = zran[1]-zran[0]
  cbform = '(E0.2)'

  mapscl = cgimgscl(mg2map_fine,minval=zran[0],maxval=zran[1],$
                    stretch=4,constant=logconstant)
  loadcv,121,/noqual,/silent
  cgimage,mapscl,/keep,pos=[0.05,0.1,0.9,0.95],opos=truepos,/noerase
  cgplot,[0],xsty=5,ysty=5,pos=truepos,/nodata,/noerase,$
         xran=[0,dx*samplefac],yran=[0,dy*samplefac]
  contourlevels = [0.02,0.04,0.08,0.16,0.32,0.64]
  cgcontour,cmap_fine,dindgen(dx*samplefac)+0.5,$
            dindgen(dy*samplefac)+0.5,$
            /overplot,c_thick=1,c_color='white',$
            levels=contourlevels,label=0
  contourlevels = [0.16,0.32,0.64]
  cgcontour,mg2map_fine,dindgen(dx*samplefac)+0.5,$
            dindgen(dy*samplefac)+0.5,$
            /overplot,c_thick=1,c_color='Black',$
            levels=contourlevels,label=0
  cbpos=[truepos[2],truepos[1],truepos[2]+0.02,truepos[3]]
  zran *= maxmg2map
  dzran *= maxmg2map
  ticknames = string(dindgen(ncbdiv+1)*dzran/double(ncbdiv) - $
                     (dzran - zran[1]),format=cbform)
   cgplot,[0],xsty=5,ysty=5,pos=truepos,/nodata,/noerase,$
          xran=[0,dx],yran=[0,dy]
  ifsf_plotaxesnuc,xran_kpc,yran_kpc,center_nuclei_kpc_x,center_nuclei_kpc_y
  cgcolorBar,position=cbpos,ticknames=ticknames,divisions=ncbdiv,$
             charsize=0.6,/right,/ver
  cgtext,0.45,0.02,'$\Delta$x (kpc)',align=0.5,/normal
  cgtext,0.05,0.55,'$\Delta$y (kpc)',align=0.5,/normal,orient=90
  cgtext,0.96,0.57,'F (MgII, erg s$\up-1$ cm$\up-2$ arcsec$\up2$)',$
          orient=270,/normal,align=0.5,chars=1

  cgps_close
  

; J1107

  zsys = 0.4665
  label = 'hizeaj1107'
  infile = dir+'cubes/'+label+'/'+label+'_3exp.fits'

  ; Planck 2018 parameters: https://ui.adsabs.harvard.edu/#abs/arXiv:1807.06209
  ldist = lumdist(zsys,H0=67.4d,Omega_m=0.315d,Lambda0=0.685d,/silent)
  asdist = ldist/(1d + zsys)^2d
  kpc_per_as = asdist*1000d/206265d
  kpc_per_pix = platescale * kpc_per_as


  cube = ifsf_readcube(infile,/quiet,header=header,datext=-1,varext=1,dqext=2)

  size_tmp = size(cube.dat)
  dx = size_tmp[1]
  dy = size_tmp[2]
  dz = size_tmp[3]
  center_axes = [double(dx)/2d,double(dy)/2d]+0.5d
  center_nuclei = [30d,30d]

;  coordinates in kpc
   xran_kpc = double([-(center_nuclei[0]-0.5),dx-(center_nuclei[0]-0.5)]) $
              * kpc_per_pix
   yran_kpc = double([-(center_nuclei[1]-0.5),dy-(center_nuclei[1]-0.5)]) $
              * kpc_per_pix
   center_nuclei_kpc_x = 0d
   center_nuclei_kpc_y = 0d

; For IDL-Colorbars to work properly
  device, true=24,decompose=0,retain=2
  ;  For interpolation + subsampling
  xresamp = (dindgen(dx*samplefac)-samplefac/2d + 0.5d)/samplefac
  yresamp = (dindgen(dy*samplefac)-samplefac/2d + 0.5d)/samplefac

  vlim = [-300,1700]
  redwave = (1d + zsys)*((linelistmg2['MgII2795']+linelistmg2['MgII2802'])/2d)
  wavesum = redwave + vlim/c_kms * redwave
  print,wavesum
  wavesub = [[wavesum[0]-100d,wavesum[0]-35d],$
             [wavesum[1]+0d,wavesum[1]+10d]]
  mg2map = ifsr_makelinemap(cube,wavesum,wavesub=wavesub,/allow)
  maxmg2map = max(mg2map.dat)
  nmg2map = mg2map.dat / maxmg2map
  mg2map_fine = interpolate(nmg2map,xresamp,yresamp,/double,/grid,cubic=-0.5)
  maxmg2map /= platescale^2d
  maxmg2map *= 1d-16

  cmap = ifsr_makelinemap(cube,[3500,5500])
  maxcmap = max(cmap.dat)
  ncmap = cmap.dat / maxcmap
  cmap_fine = interpolate(ncmap,xresamp,yresamp,/double,/grid,cubic=-0.5)
  maxcmap /= (0.2914d)^2d
  maxcmap *= 1d-16


  cgps_open,dir+'maps/'+label+'/rb1/'+label+'mg2.eps',charsize=1d,$
            /encap,/inches,xs=7.5d,ys=7.5d,/qui,/nomatch

  zran = [0.1,1]
  ncbdiv = 5
  dzran = zran[1]-zran[0]
  cbform = '(E0.2)'

  mapscl = cgimgscl(mg2map_fine,minval=zran[0],maxval=zran[1],$
                    stretch=4,constant=logconstant)
  loadcv,121,/noqual,/silent
  cgimage,mapscl,/keep,pos=[0.05,0.1,0.9,0.95],opos=truepos,/noerase
  cgplot,[0],xsty=5,ysty=5,pos=truepos,/nodata,/noerase,$
         xran=[0,dx*samplefac],yran=[0,dy*samplefac]
  contourlevels = [0.02,0.04,0.08,0.16,0.32,0.64]
  cgcontour,cmap_fine,dindgen(dx*samplefac)+0.5,$
            dindgen(dy*samplefac)+0.5,$
            /overplot,c_thick=1,c_color='white',$
            levels=contourlevels,label=0
  contourlevels = [0.16,0.32,0.64]
  cgcontour,mg2map_fine,dindgen(dx*samplefac)+0.5,$
            dindgen(dy*samplefac)+0.5,$
            /overplot,c_thick=1,c_color='Black',$
            levels=contourlevels,label=0
  cbpos=[truepos[2],truepos[1],truepos[2]+0.02,truepos[3]]
  zran *= maxmg2map
  dzran *= maxmg2map
  ticknames = string(dindgen(ncbdiv+1)*dzran/double(ncbdiv) - $
                     (dzran - zran[1]),format=cbform)
   cgplot,[0],xsty=5,ysty=5,pos=truepos,/nodata,/noerase,$
          xran=[0,dx],yran=[0,dy]
  ifsf_plotaxesnuc,xran_kpc,yran_kpc,center_nuclei_kpc_x,center_nuclei_kpc_y
  cgcolorBar,position=cbpos,ticknames=ticknames,divisions=ncbdiv,$
             charsize=0.6,/right,/ver
  cgtext,0.45,0.02,'$\Delta$x (kpc)',align=0.5,/normal
  cgtext,0.05,0.55,'$\Delta$y (kpc)',align=0.5,/normal,orient=90
  cgtext,0.96,0.57,'F (MgII, erg s$\up-1$ cm$\up-2$ arcsec$\up2$)',$
          orient=270,/normal,align=0.5,chars=1

  cgps_close
  
end
