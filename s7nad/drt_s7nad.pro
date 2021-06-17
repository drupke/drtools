; Have to .COMPILE PPXF beforehand to access helper routine PPXF_REDDENING_CURVE
; RESOLVE_ROUTINE statements below fail, though PPXF_REDDENING_CURVE is 
; recognized on subsequent runs of DRT_S7NAD. Weird.

pro drt_s7nad,redolinmix=redolinmix,redofitexy=redofitexy

;  So that IDL can find the helpder routine PPXF_REDDENING_CURVE
;   resolve_routine, 'ppxf'
;   resolve_routine, 'ppxf_reddening_curve',/is_function

   gals = ['ngc1266','ngc1808','eso500','ngc5728',$
           'eso339','ic5063','ic5169','ic1481']

   plotdir = '/Users/drupke/ifs/s7/plots/s7nad/'
   mapdir = '/Users/drupke/ifs/s7/maps/'
   specfitdir = '/Users/drupke/specfits/s7/'
   tabdir = '/Users/drupke/Documents/papers_and_talks/papers/s7/tabs/'

   amp = '&'
   dslash = '\\'
   ndat = '\nodata'
   lineofdashes = strjoin(replicate('-',62))

   substar = '!D!10'+string(72B)+'!X!N'

   plotquantum = 2.5d ; in inches
   bad = 1d99
   
;  factor by which to resample images for PS-to-PDF conversion
   samplefac = 10
   samplefac2 = 100

;  # of bootstrap samples for FITEXY errors
   nsamp = 1000
   lm_miniter = 7894L; 15787 is 4sigma, and miniter typically gets doubled

;  arrays to hold disk areas, fractional inflow/outflow areas, and sizes
   maxrad_disk = dblarr(8)
   maxrad_in = dblarr(8)
   maxrad_out = dblarr(8)
   diskarea_nspax = intarr(8)
   diskarea_sqkpc = dblarr(8)
   infracarea = dblarr(8)
   outfracarea = dblarr(8)
   
;  arrays to hold averages
   nnaiall_in = []
   nnaiall_out = []
   nnaiavg_in = dblarr(8)
   nnaiavg_out = dblarr(8)
   nnaisdev_in = dblarr(8)
   nnaisdev_out = dblarr(8)
   cfall_in = []
   cfall_out = []
   cfavg_in = dblarr(8)
   cfavg_out = dblarr(8)
   cfsdev_in = dblarr(8)
   cfsdev_out = dblarr(8)
   weqall_in = []
   weqall_out = []
   weqavg_in = dblarr(8)
   weqavg_out = dblarr(8)
   weqsdev_in = dblarr(8)
   weqsdev_out = dblarr(8)


   tvlct,[[27],[158],[119]],100
   tvlct,[[217],[95],[2]],101
   tvlct,[[117],[112],[179]],102
   tvlct,[[231],[41],[138]],103
   tvlct,[[0],[0],[0]],104
   tvlct,[[102],[166],[30]],105

;  areas and Weq(abs), Cf/tau corr. coeffs, stellar errors
   sqkpc_per_spaxel = dblarr(n_elements(gals))
   weq_abs_cf_cc = dblarr(n_elements(gals))
   weq_abs_tau_cc = dblarr(n_elements(gals))
   med_stel_errvel = dblarr(n_elements(gals))
   med_stel_errsig = dblarr(n_elements(gals))
   med_stel_errebv = dblarr(n_elements(gals))
   cf_stats = dblarr(n_elements(gals),3)
   weq_abs_stats = dblarr(n_elements(gals),3)
   for i=0,n_elements(gals)-1 do begin
      xdr = mapdir+gals[i]+'/vormerge/'+gals[i]+'_plots.xdr'
      if file_test(xdr) then begin
         restore,file=xdr
         kpc_per_as = plotstr.kpc_per_as
         sqkpc_per_spaxel[i] = kpc_per_as^2d
         weq_abs_cf_cc[i] = plotstr.weq_abs_cf_cc
         weq_abs_tau_cc[i] = plotstr.weq_abs_tau_cc
         igdstel = where(plotstr.stel_vel ne bad)
         stel_errvel = $
            (plotstr.stel_errvel[*,*,0] + plotstr.stel_errvel[*,*,1])/2d
         stel_errsig = $
            (plotstr.stel_errsig[*,*,0] + plotstr.stel_errsig[*,*,1])/2d
         stel_errebv = $
            (plotstr.stel_errebv[*,*,0] + plotstr.stel_errebv[*,*,1])/2d
         ; get rid of duplicates due to Voronoi tiling
         errvelall = stel_errvel[igdstel]
         errsigall = stel_errsig[igdstel]
         errebvall = stel_errebv[igdstel]
         zall = errvelall+errsigall+errebvall
         isort = bsort(zall)
         iuniq = uniq(zall[isort])
         isortuniq = isort[iuniq]
         errvel_unsort = errvelall[isortuniq]
         errsig_unsort = errsigall[isortuniq]
         errebv_unsort = errebvall[isortuniq]
         iresort = bsort(errvel_unsort)
         stel_errvel = errvel_unsort[iresort]
         stel_errsig = errsig_unsort[iresort]
         stel_errebv = errebv_unsort[iresort]
         ; ... end getting rid of duplicates
         med_stel_errvel[i] = median(stel_errvel)
         med_stel_errsig[i] = median(stel_errsig)
         med_stel_errebv[i] = median(stel_errebv)
         ; cf stats
         cf_stats[i,*] = plotstr.cf_stats
         weq_abs_stats[i,*] = plotstr.weq_abs_stats
      endif
   endfor
   ingc1808 = where(gals eq 'ngc1808')
   ref_sqkpc_per_spaxel = sqkpc_per_spaxel[where(gals eq 'ngc1808')]
   areanorm = sqkpc_per_spaxel/ref_sqkpc_per_spaxel[0]
   intareanorm = round(areanorm)

   ; array of all sigmas
   allsig_abs = []
   allsig_em = []
   for i=0,7 do begin
      xdr = mapdir+gals[i]+'/vormerge/'+gals[i]+'_plots.xdr'
      if file_test(xdr) then begin
         restore,file=xdr
         initdat = call_function('ifsf_'+gals[i]+'vormerge')
         vormap=initdat.vormap
         nvorcols = max(vormap)
         vorcoords = intarr(nvorcols,2)
         for j=1,nvorcols do begin
            ivor = where(vormap eq j,ctivor)
            xyvor = array_indices(vormap,ivor[0])
            size_abs = size(plotstr.sig_abs)
            if size_abs[0] gt 1 then $
               if plotstr.sig_abs[xyvor[0],xyvor[1]] ne bad then $
                  allsig_abs = [allsig_abs, plotstr.sig_abs[xyvor[0],xyvor[1]]]
            size_em = size(plotstr.sig_em)
            if size_em[0] gt 1 then $
               if plotstr.sig_em[xyvor[0],xyvor[1]] ne bad then $
                  allsig_em = [allsig_em, plotstr.sig_em[xyvor[0],xyvor[1]]]

         endfor
      endif
   endfor


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; NaD: N(NaI) and Weq_total vs. E(B-V): other data
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

   zeroebv = -99d

;  M31
   ebv_m31 = [0.07,0.32,0.33,0.00,0.24,0.05,0.09,-0.03,0.28,0.22,$
              0.17,-0.01,0.05,0.18,0.06,0.43,0.13]
   ebv_err_m31 = [0.04,0.04,0.04,0.04,0.04,0.04,0.04,0.04,0.04,0.05,$
                  0.04,0.04,0.04,0.06,0.05,0.04,0.05]
   nnai_m31 = [12.90,13.57,13.16,12.36,13.29,13.56,12.24,11.79,12.78,13.08,$
               13.26,11.4,13.48,13.10,13.16,13.25,13.22]
   nnai_errhi_m31 = [0.14,0.38,0.11,0.56,0.45,0.63,0.03,0.12,0.02,0.22,$
                     0.51,0.00,0.62,0.95,0.18,0.05,0.47]
   nnai_errlo_m31 = [0.13,0.29,0.10,0.33,0.41,0.42,0.03,0.10,0.02,0.35,$
                     0.35,11.4,0.24,0.36,0.10,0.03,0.27]

   logebv_m31 = alog10(ebv_m31)   
   logebv_errlo_m31 = logebv_m31 - alog10(ebv_m31 - ebv_err_m31)
   logebv_errhi_m31 = alog10(ebv_m31 + ebv_err_m31) - logebv_m31
   
   izeroebv = where(ebv_m31 le 0)
   logebv_m31[izeroebv] = zeroebv

;  SNe
;  upper limits not included.
   nnai_sne = [13.740,11.800,12.970,12.512,12.76,14.406,12.329,12.886,13.246,$
               13.779,15.242,12.044,12.844,12.583,12.595,13.281,13.220,12.720,$
               13.055,14.472,12.920,13.105,13.254,12.605,12.701,11.495,12.989]
   nnai_err_sne = [.050,.020,.470,.073,.03,.862,.088,.125,.021,.041,.069,.142,$
                   .024,.029,.171,.012,.182,.184,.044,.036,.058,.161,.103,$
                   .021,.029,.024,.035]
   rv_sne = [2.57,2.24,2.28,1.82,2.25,1.22,2.16,2.45,2.46,1.31,1.95,2.11,2.17,$
            2.24,2.13,1.46,1.74,2.42,2.21,1.20,2.20,2.29,2.48,2.24,1.54,1.63,2.7]
   rv_errhi_sne = [.23,.62,.62,.76,.46,.26,.66,.52,.57,.08,.46,.55,.53,.64,$
                   .55,.32,.50,.56,.53,.26,.55,.68,.57,.75,.57,.60,.5]
   rv_errlo_sne = [.21,.73,.80,.53,.36,.21,.68,.77,.63,.10,.33,.48,.60,.78,$
                   .64,.24,.64,.72,.73,.14,.69,.62,.83,.74,.59,.53,.5]
   av_sne = [2.03,.08,.17,.49,.62,.62,.21,.08,.17,1.88,2.15,.31,.09,.12,$
             .04,.54,.23,.39,.53,.71,.16,.05,.24,.39,.50,.15,.67]
   av_errhi_sne = [.09,0,0,.19,.08,.09,.09,0,0,.09,.40,.06,0,.08,$
                   0,.08,.06,.16,.13,.10,.08,0,.09,.14,.17,.19,.12]
   av_errlo_sne = [.13,.08,.17,.14,.08,.10,.07,.08,.17,.13,.35,.06,.09,.05,$
                   .04,.08,.07,.11,.16,.08,.05,.05,.08,.11,.19,.06,.12]
;  Av upper limits have no EW err listed, so assume 0.1. No EW listed for 1994D.
   weqnai_sne = [3.91,0.00,0.59,0.56,0.67,2.75,0.35,0.90,2.25,$
                 2.35,2.23,0.31,0.79,0.45,0.74,1.55,0.35,0.99,$
                 0.76,2.32,1.20,0.50,1.95,0.72,0.46,0.07,1.59]
   weqnai_err_sne = [0.13,0.00,0.1,0.14,0.08,0.10,0.07,0.1,0.27,$
                     0.13,0.35,0.06,0.1,0.05,0.1,0.08,0.07,0.11,$
                     0.16,0.08,0.05,0.1,0.08,0.11,0.19,0.06,0.12]
   n_sne = n_elements(nnai_sne)

   ebv_sne = av_sne/rv_sne
   ebv_errlo_sne = ebv_sne
   ebv_errhi_sne = dblarr(n_sne)
   ebv_det_sne = where(av_errhi_sne ne 0)
   ebv_errlo_sne[ebv_det_sne] = $
      ebv_sne[ebv_det_sne]*$
      sqrt((av_errlo_sne[ebv_det_sne]/av_sne[ebv_det_sne])^2d + $
           (rv_errhi_sne[ebv_det_sne]/rv_sne[ebv_det_sne])^2d)
   ebv_errhi_sne[ebv_det_sne] = $
      ebv_sne[ebv_det_sne]*$
      sqrt((av_errhi_sne[ebv_det_sne]/av_sne[ebv_det_sne])^2d + $
           (rv_errlo_sne[ebv_det_sne]/rv_sne[ebv_det_sne])^2d)

   logebv_sne = alog10(ebv_sne)
   logebv_errlo_sne = dblarr(n_sne) + 99d
   logebv_errhi_sne = dblarr(n_sne)
   logebv_errlo_sne[ebv_det_sne] = $
      logebv_sne[ebv_det_sne] - $
      alog10(ebv_sne[ebv_det_sne] - ebv_errlo_sne[ebv_det_sne])
   logebv_errhi_sne[ebv_det_sne] = $
      alog10(ebv_sne[ebv_det_sne] + ebv_errhi_sne[ebv_det_sne]) - $
      logebv_sne[ebv_det_sne]

   logweqnai_sne = alog10(weqnai_sne)
   logweqnai_errlo_sne = dblarr(n_sne) + 99d
   logweqnai_errhi_sne = dblarr(n_sne)
   weqnai_det_sne = where(weqnai_sne ne 0)
   logweqnai_errlo_sne[weqnai_det_sne] = $
      logweqnai_sne[weqnai_det_sne] - $
      alog10(weqnai_sne[weqnai_det_sne] - weqnai_err_sne[weqnai_det_sne])
   logweqnai_errhi_sne[weqnai_det_sne] = $
      alog10(weqnai_sne[weqnai_det_sne] + weqnai_err_sne[weqnai_det_sne]) - $
      logweqnai_sne[weqnai_det_sne]
           
;  LMC (Cox + 06)
   nnai_lmc = [12.59,12.84,13.49,13.25,12.49,11.30]
   nnai_err_lmc = [0.03,0.03,0.08,0.08,0.03,0.03]        
   ebv_lmc = [.20,.29,.12,.14,.04,.02]
   ebv_err_lmc = [.06,.04,.04,.02,.03,.02]
   
   logebv_lmc = alog10(ebv_lmc)
   logebv_errlo_lmc = logebv_lmc - alog10(ebv_lmc - ebv_err_lmc)
   logebv_errhi_lmc = alog10(ebv_lmc + ebv_err_lmc) - logebv_lmc
          
;  LMC/SMC (Welty + 06)
   nnai_mc = [12.81,12.98,11.96,12.15,13.72,13.38,11.92,12.86,13.79,12.50,13.0,$
              12.69,13.37,13.17,12.42,12.59,13.0,12.55,12.84,13.01,12.68]
   ebv_mc = [.14,.08,.07,.12,.34,.18,.11,.26,.22,.13,.12,.09,.51,.21,.10,.20,$
             .13,.36,.38,.08,.14]

   logebv_mc = alog10(ebv_mc)
             
;  U/LIRGs
;   F1_5,F09320+6134,F10190+1322:W,F10190+1322:E,F15386+3807,F16333+4630:E,
;   F16474+3430:N,most ULIRGs from Rupke+05 ...
   nnai_lirgs = [13.73,13.77,13.89,13.92,14.31,13.03,14.29,13.35,14.76,13.41,$
                 13.34,13.45,13.20,13.30,14.17,13.63,13.58,13.19,14.25,13.67,$
                 13.78,13.89,14.01,13.68,13.61,13.44]
   ebv_lirgs = [0.91,1.63,0.89,1.88,1.21,0.49,0.57,0.92,0.90,1.66,0.75,$
                0.79,1.56,0.69,0.82,0.54,1.06,0.97,1.09,1.43,0,1.02,1.13,$
                0.58,0.82,0.35,0.65]
   weq_lirgs = [5.18,3.55,5.67,8.85,5.99,1.81,2.62,4.82,6.90,2.42,3.20,2.76,$
                1.99,5.15,3.33,3.97,5.05,3.42,4.05,2.43,7.85,2.66,2.13,4.24,$
                2.28,4.87]
                
   logebv_lirgs = alog10(ebv_lirgs)
   
   cgps_open,plotdir+'s7nad_nnai_v_ebv_gas_other.eps',charsize=1,$
             /encap,/nomatch,/inches,xs=7.5d,ys=7.5d,/qui
;  units are pixels per tenths of an inch
   xran = [-2,0.5]
   yran = [11,15.5]
   cgplot,[0],/xsty,/ysty,xran=xran,yran=yran,pos=[0.2,0.1,0.95,0.95],$
          ytickform='(A1)'

; Fitted lines

;   ;  Model of screen absorption (?) and constant gas-to-dust ratio
;   av=dindgen(21)/20d*1.75d ; gives max. N(NaI) of 13.06, which correponds to
;   ; max of N(H) 21.5 observed in Wakker & Mathis
;   ;  NH/AV from Predehl & Schmitt 1995
;   nh=alog10(av) + 21d + alog10(1.8d)
;   ;  Wakker & Mathis 2000:
;   ;  log [N(NaI)/N(HI)] = (-0.16+/-0.06) log N(HI) - (5.00+/-4.94)
;   nnai = nh - 0.16d*nh - 5d
;   logebv_rv325 = alog10(av/3.25d)
;   logebv_rv405 = alog10(av/4.05d)
;   logebv_rv485 = alog10(av/4.85d)
;   cgoplot,logebv_rv325,nnai,linesty=2
;   cgoplot,logebv_rv405,nnai,linesty=2
;   cgoplot,logebv_rv485,nnai,linesty=2

;  Hobbs 1974
;  N(NaI) = 1.7x10^14 E(B-V)^1.8
   ebv_hobbs = dindgen(26)/25d*0.5d
   logebv_hobbs = alog10(ebv_hobbs)
   nnai_hobbs = 14d + alog10(1.7) + 1.8d*logebv_hobbs
   nnai_hobbs_extrap = 14d + alog10(1.7) + 1.8d*(logebv_hobbs+1d)
   cgoplot,logebv_hobbs,nnai_hobbs,color=105,thick=4
   cgoplot,logebv_hobbs,nnai_hobbs+0.3d,color=105,thick=2
   cgoplot,logebv_hobbs,nnai_hobbs-0.3d,color=105,thick=2
   cgoplot,logebv_hobbs+1d,nnai_hobbs_extrap,color=105,thick=4,linesty=2
   cgoplot,logebv_hobbs+1d,nnai_hobbs_extrap+0.3d,color=105,thick=2,linesty=2
   cgoplot,logebv_hobbs+1d,nnai_hobbs_extrap-0.3d,color=105,thick=2,linesty=2

;  Phillips et al. 2013: fit to Sembach+93 + their own MW data
   ebv_sembach_phillips = dindgen(41)/40d*0.4d
   logebv_sembach_phillips = alog10(ebv_sembach_phillips)
   nnai_sembach_phillips = 13.180d + 1.125d*(logebv_sembach_phillips+alog10(3d))
   nnai_sembach_phillips_extrap = 13.180d + 1.125d*(logebv_sembach_phillips+alog10(3d)+1d)
   cgoplot, logebv_sembach_phillips, nnai_sembach_phillips,color='Black',thick=4
   cgoplot, logebv_sembach_phillips, nnai_sembach_phillips+0.26d,color='Black',thick=2
   cgoplot, logebv_sembach_phillips, nnai_sembach_phillips-0.26d,color='Black',thick=2
   cgoplot, logebv_sembach_phillips+1d, nnai_sembach_phillips_extrap,color='Black',thick=4,linesty=2
   cgoplot, logebv_sembach_phillips+1d, nnai_sembach_phillips_extrap+0.26d,color='Black',thick=2,linesty=2
   cgoplot, logebv_sembach_phillips+1d, nnai_sembach_phillips_extrap-0.26d,color='Black',thick=2,linesty=2

;  Baron et al. 2016:
;  N(NaI) = (3.4+/-0.2)*10^13 E(B-V)
;  Don't include because optically thin regime assumed; no measurement of N(NaI)
;  Line matches Phillips+2013 fit
;   ebv_baron = dindgen(31)/30d*0.3d
;   logebv_baron = alog10(ebv_baron)
;   nnai_baron = 13d + alog10(3.4)+logebv_baron
;   cgoplot,logebv_baron,nnai_baron,linesty=3,color=103,thick=8
;   cgoplot,logebv_baron,nnai_baron+0.3,linesty=3,color=103,thick=8
;   cgoplot,logebv_baron,nnai_baron-0.3,linesty=3,color=103,thick=8

   xtit='log [ Gas E(B-V) ]'
   ytit='log [ N(NaI) / cm$\up-2$ ]'
   cgtext,xtit,0.5d,0.02d,align=0.5,/norm,chars=1.5
   cgtext,ytit,0.03d,0.5d,align=0.5,orient=90d,/norm,chars=1.5

;errors
;  M31 pionts
   cgoplot,logebv_m31,nnai_m31,err_xlow=logebv_errlo_m31,err_xhi=logebv_errhi_m31,$
           err_ylow=nnai_errlo_m31,err_yhi=nnai_errhi_m31,/err_clip,$
           psym=1,err_color='Grey',err_width=0
;  MC pionts
   cgoplot,logebv_lmc,nnai_lmc,err_xlow=logebv_errlo_lmc,err_xhi=logebv_errhi_lmc,$
           err_ylow=nnai_errlo_lmc,err_yhi=nnai_errhi_lmc,/err_clip,$
           psym=1,err_color='Grey',err_width=0
;  SNe points
   cgoplot,logebv_sne,nnai_sne,err_xlow=logebv_errlo_sne,err_xhi=logebv_errhi_sne,$
           err_ylow=nnai_err_sne,err_yhi=nnai_err_sne,/err_clip,$
           psym=1,err_color='Grey',err_width=0
;data
;  M31 pionts
   cgoplot,logebv_m31,nnai_m31,psym=15,color=100,symsize=1
   cgoplot,logebv_m31,nnai_m31,psym=6,symsize=1
;  MC pionts
   cgoplot,logebv_lmc,nnai_lmc,psym=15,color=101,symsize=1
   cgoplot,logebv_lmc,nnai_lmc,psym=6,symsize=1
   cgoplot,alog10(ebv_mc),nnai_mc,psym=15,color=101,symsize=1
   cgoplot,alog10(ebv_mc),nnai_mc,psym=6,symsize=1
;  SNe points
   cgoplot,logebv_sne,nnai_sne,psym=15,color=102,symsize=1
   cgoplot,logebv_sne,nnai_sne,psym=6,symsize=1
;  LIRGs
   cgoplot,logebv_lirgs,nnai_lirgs,psym=15,color=103,symsize=1
   cgoplot,logebv_lirgs,nnai_lirgs,psym=6,symsize=1

   al_legend,['M31 (Cordiner+2011)','LMC/SMC (Cox+2006, Welty+2006)',$
              'SNe (Phillips+2013)','U/LIRGs (Rupke+2005,2008)'],$
              color=[100,101,102,103],psym=[15,15,15,15],$
             symsize=[1,1,1,1],/top,/left
   al_legend,['MW (Hobbs 1974)','MW (Phillips+2013)'],$ ;,'QSOs (Baron+2016)'],$
;             color=[105,104,103],linesty=[2,0,3],/bottom,/right,$
;             thick=[8,4,8]
             color=[105,104],linesty=[0,0],/bottom,/right,$
             thick=[4,4]

   cgoplot,[-1.7,0.3],[13.05,15.55],color=101,thick=8,linesty=1
   cgoplot,[-1.7,0.3],[10.6,13.35],color=101,thick=8,linesty=1


;  E(B-V) = 0
   xran = [-99.1,-98.9]
   yran = [11,15.5]
   cgplot,[0],/xsty,/ysty,xran=xran,yran=yran,pos=[0.1,0.1,0.2,0.95],/noerase,$
          xticks=2,xtickn=[' ',' ',' '],yticklen=!P.ticklen*10
   cgoplot,logebv_m31,nnai_m31,err_xlow=logebv_errlo_m31,err_xhi=logebv_errhi_m31,$
           err_ylow=nnai_errlo_m31,err_yhi=nnai_errhi_m31,/err_clip,$
           psym=1,err_color='Grey',err_width=0
   cgoplot,logebv_m31,nnai_m31,psym=15,color=100,symsize=1
   cgoplot,logebv_m31,nnai_m31,psym=6,symsize=1

           
   cgps_close
   
   cgps_open,plotdir+'s7nad_weq_v_ebv_gas_other.eps',charsize=1,$
             /encap,/nomatch,/inches,xs=7.5d,ys=7.5d,/qui
;  units are pixels per tenths of an inch
   xran = [-1.8,0.3]
   yran = [-1,1]
   xtit='log [ Gas E(B-V) ]'
   ytit='log [ W$\downeq$ / $\angstrom$ ]'
   cgplot,[0],/xsty,/ysty,xran=xran,yran=yran,pos=[0.15,0.1,0.95,0.95],$
          xtit=xtit,ytit=ytit,charsize=1.5

;;  from Richmond et al. 1994
;  A lot of scatter, don't use ...
;   ebv_modr = dindgen(11)/10d*0.4d
;   weq_modr = -0.0956d + 4.91d*ebv_modr
;   cgoplot,alog10(ebv_modr),alog10(weq_modr),color=101
;  relationship from Poznanski et al. 2012; scatter is negligible (80 mA)
   ebv_modp = dindgen(101)/100d*2d
   logebv_modp = alog10(ebv_modp)
   weq_modp = 1/1.17d*(logebv_modp+1.85d)
   cgoplot,logebv_modp,alog10(weq_modp),color='Black',thick=4
   cgoplot,logebv_modp,alog10(weq_modp+0.080d),color='Black',thick=2,linesty=2
   cgoplot,logebv_modp,alog10(weq_modp-0.080d),color='Black',thick=2,linesty=2
;  relationship from Baron et al. 2016
   ebv_modb = dindgen(11)/10d*0.3d
   weq_modb = 10.08d*ebv_modb
   cgoplot,alog10(ebv_modb),alog10(weq_modb),color=105,thick=4
   cgoplot,alog10(ebv_modb),alog10(weq_modb+0.55d),color=105,thick=2,linesty=2
   cgoplot,alog10(ebv_modb),alog10(weq_modb-0.55d),color=105,thick=2,linesty=2
;  relationship from Chen et al. 2010, estimated
   ebv_modb = dindgen(11)/10d*(2.5d - 1.3d)/3.1d + 1.3d/3.1d
   weq_modb = 1d/(1.3d/3.1d)*ebv_modb
   cgoplot,alog10(ebv_modb),alog10(weq_modb),color=102,thick=4

;  SNe points
   cgoplot,logebv_sne,logweqnai_sne,err_xlow=logebv_errlo_sne,err_xhi=logebv_errhi_sne,$
           err_ylow=logweqnai_errlo_sne,err_yhi=logweqnai_errhi_sne,/err_clip,$
           psym=1,err_color='Grey',err_width=0
;  SNe points
   cgoplot,logebv_sne,logweqnai_sne,psym=15,color=102,symsize=1
   cgoplot,logebv_sne,logweqnai_sne,psym=6,symsize=1
;  LIRGs
   cgoplot,logebv_lirgs,alog10(weq_lirgs),psym=15,color=103,symsize=1
   cgoplot,logebv_lirgs,alog10(weq_lirgs),psym=6,symsize=1

   cgoplot,[-1.7,0.3],[0.1,1.05],color=101,thick=8,linesty=1
   cgoplot,[-1.7,0.3],[-1.4,0.3],color=101,thick=8,linesty=1


   al_legend,['SNe (Phillips+2013)','U/LIRGs (Rupke+2005,2008)'],$
              color=[102,103],psym=[15,15],$
             symsize=[1,1],/top,/left
   al_legend,['Galaxies (Chen+2010)',$
              'MW (Poznanski+2012)','AGN (Baron+2016)'],$ ;,'QSOs (Baron+2016)'],$
             color=[102,104,105],linesty=[0,0,0],/bottom,/right,$
             thick=[4,4,4]
           
   cgps_close
   
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; NaD: Total Weq vs. Stellar E(B-V)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
     
   nx = 4
   ny = 2
   ngal = nx*ny

   if keyword_set(redolinmix) then redolinmix_totweq_v_ebv_stel = 1b $
   else redolinmix_totweq_v_ebv_stel = 0b
   ; redolinmix_totweq_v_ebv_stel = 1b
   if keyword_set(redofitexy) then redofitexy_totweq_v_ebv_stel = 1b $
   else redofitexy_totweq_v_ebv_stel = 0b
   ; redofitexy_totweq_v_ebv_stel = 1b

   openw,lun_tmp,tabdir+'table_totweq_v_ebv_stel.txt',/get_lun
 
   cgps_open,plotdir+'s7nad_totweq_v_ebv_stel.eps',charsize=1,$
             /encap,/nomatch,/inches,xs=plotquantum*nx*0.95d,ys=plotquantum*ny,/qui
;  units are pixels per tenths of an inch
   pos = cglayout([nx,ny],oxmar=[6,0],oymar=[6,3],ixmar=[-1,0],$
                  unit=!D.X_PX_CM*2.54d/10d,xgap=4d,ygap=5d)
   noerase=0b
   fit1_rxyall=dblarr(ngal)
   fit1_weq0all=dblarr(ngal)
   fit1_mall_linmix=dblarr(ngal)
   fit1_mall_fitexy=dblarr(ngal)
   fit1_dweqall_linmix=dblarr(ngal)
   fit1_dweqall_fitexy=dblarr(ngal)
   maxpts = 1000 ; maximum number of Voro
   if redolinmix_totweq_v_ebv_stel then begin
      linmixpar = {alpha:dblarr(ngal,5),$
                   beta:dblarr(ngal,5),$
                   corr:dblarr(ngal,5),$
                   sigsqr:dblarr(ngal,5),$
                   nfit:dblarr(ngal),$
                   pval:dblarr(ngal,2),$
                   xdat:dblarr(ngal,maxpts)+bad,$
                   ydat:dblarr(ngal,maxpts)+bad,$
                   xerr:dblarr(ngal,maxpts)+bad,$
                   yerr:dblarr(ngal,maxpts)+bad,$
                   yfit:dblarr(ngal,maxpts)+bad,$
                   yfitlo:dblarr(ngal,maxpts)+bad,$
                   yfithi:dblarr(ngal,maxpts)+bad,$
                   yfit2lo:dblarr(ngal,maxpts)+bad,$
                   yfit2hi:dblarr(ngal,maxpts)+bad}
   endif else begin
      fitfile = plotdir+'totweq_v_ebv_stel_linmix.xdr'
      if file_test(fitfile) then $
         restore,file=fitfile $
      else stop,'WARNING: File '+fitfile+' does not exist; aborting.'
   endelse
   if redofitexy_totweq_v_ebv_stel then begin
      fitexypar = {result:dblarr(ngal,2),$
         nfit:intarr(ngal),$
         ymod:dblarr(ngal,maxpts)+bad,$
         scat:dblarr(ngal),$
         berr:dblarr(ngal),$
         merr:dblarr(ngal),$
         scaterr:dblarr(ngal)}
   endif else begin
      fitfile = plotdir+'totweq_v_ebv_stel_fitexy.xdr'
      if file_test(fitfile) then $
         restore,file=fitfile $
      else stop,'WARNING: File '+fitfile+' does not exist; aborting.'
   endelse
   for i=0,ngal-1 do begin
      xdr = mapdir+gals[i]+'/vormerge/'+gals[i]+'_plots.xdr'
      if file_test(xdr) then begin
         restore,file=xdr
         x = plotstr.stel_ebv
         xerrlo = plotstr.stel_errebv[*,*,0]
         xerrhi = plotstr.stel_errebv[*,*,1]
         y = plotstr.weq_abs_plus_em
         yboth = plotstr.weq_abs_and_em
         yerrlo = plotstr.weqerrlo_abs_plus_em
         yerrhi = plotstr.weqerrhi_abs_plus_em
         igd = where(x ne bad AND y ne bad)
;         xran = [0.95d,1.05d]*[min(x[igd]),max(x[igd])]
;         if xran[0] le 0d then xran[0] = -0.05d*max(x[igd])
;         yran = [min(y[igd]),max(y[igd])]
;         if yran[0] lt 0 then yran *= [1.05d,1.05d] $
;         else if yran[0] eq 0 then yran = [-0.05d*yran[1],1.05*yran[1]] $
;         else yran *= [0.95d,1.05d]
         xran = [-0.049,1.199]
         yran = [-4,11]
         if i gt 0 then noerase=1b
         if ~ redolinmix_totweq_v_ebv_stel then begin
            sigsqr = reform(linmixpar.sigsqr[i,*],5)
            alpha = reform(linmixpar.alpha[i,*],5)
            beta = reform(linmixpar.beta[i,*],5)
            corr = reform(linmixpar.corr[i,*],5)
            xdat = reform(linmixpar.xdat[i,0:linmixpar.nfit[i]-1],$
                          linmixpar.nfit[i])
            ydat = reform(linmixpar.ydat[i,0:linmixpar.nfit[i]-1],$
                          linmixpar.nfit[i])
            xerr = reform(linmixpar.xerr[i,0:linmixpar.nfit[i]-1],$
                          linmixpar.nfit[i])
            yerr = reform(linmixpar.yerr[i,0:linmixpar.nfit[i]-1],$
                          linmixpar.nfit[i])
            yfit = reform(linmixpar.yfit[i,0:linmixpar.nfit[i]-1],$
                          linmixpar.nfit[i])
            yfitlo = reform(linmixpar.yfitlo[i,0:linmixpar.nfit[i]-1],$
                            linmixpar.nfit[i])
            yfithi = reform(linmixpar.yfithi[i,0:linmixpar.nfit[i]-1],$
                            linmixpar.nfit[i])
            yfit2lo = reform(linmixpar.yfit2lo[i,0:linmixpar.nfit[i]-1],$
                            linmixpar.nfit[i])
            yfit2hi = reform(linmixpar.yfit2hi[i,0:linmixpar.nfit[i]-1],$
                            linmixpar.nfit[i])
         endif else begin
            print,'Running LINMIX for ',gals[i]
            ifit = where(x ne bad AND y ne bad)
;           get rid of duplicates due to Voronoi tiling
            xdatall = x[ifit]
            xerrall = (xerrlo[ifit]+xerrhi[ifit])/2d
            ydatall = y[ifit]
            yerrall = (yerrlo[ifit]+yerrhi[ifit])/2d
            zdatall = xdatall+ydatall
;           Unique elements from x+y
            isort = bsort(zdatall)           
            iuniq = uniq(zdatall[isort])
            isortuniq = isort[iuniq]
            xdat_unsort = xdatall[isortuniq]
            ydat_unsort = ydatall[isortuniq]
            xerr_unsort = xerrall[isortuniq]
            yerr_unsort = yerrall[isortuniq]
            iresort = bsort(xdat_unsort)
            xdat = xdat_unsort[iresort]
            ydat = ydat_unsort[iresort]
            xerr = xerr_unsort[iresort]
            yerr = yerr_unsort[iresort]
            fitout = drt_runlinmix(xdat,ydat,xerr=xerr,yerr=yerr,/metro,$
                                   alpha=alpha,beta=beta,corr=corr,$
                                   sigsqr=sigsqr,pval=pval,miniter=lm_miniter)
            print,pval
            yfit = fitout[*,0]
            yfitlo = fitout[*,1]
            yfithi = fitout[*,2]
            yfit2lo = fitout[*,3]
            yfit2hi = fitout[*,4]
            nfit = n_elements(xdat)
            linmixpar.alpha[i,*] = alpha
            linmixpar.beta[i,*] = beta
            linmixpar.corr[i,*] = corr
            linmixpar.sigsqr[i,*] = sigsqr
            linmixpar.nfit[i] = nfit
            linmixpar.pval[i,*] = pval
            linmixpar.xdat[i,0:nfit-1] = xdat
            linmixpar.ydat[i,0:nfit-1] = ydat
            linmixpar.xerr[i,0:nfit-1] = xerr
            linmixpar.yerr[i,0:nfit-1] = yerr
            linmixpar.yfit[i,0:nfit-1] = yfit
            linmixpar.yfitlo[i,0:nfit-1] = yfitlo
            linmixpar.yfithi[i,0:nfit-1] = yfithi
            linmixpar.yfit2lo[i,0:nfit-1] = yfit2lo
            linmixpar.yfit2hi[i,0:nfit-1] = yfit2hi
         endelse
         ; FITEXY fits
         if ~ redofitexy_totweq_v_ebv_stel then begin
            fitexy_result = reform(fitexypar.result[i,*],2)
            nfit = fitexypar.nfit[i]
            ymod = reform(fitexypar.ymod[i,0:nfit-1],nfit)
            fitexy_scat = fitexypar.scat[i]
            fitexy_merr = fitexypar.merr[i]
            fitexy_berr = fitexypar.berr[i]
            fitexy_scaterr = fitexypar.scaterr[i]
         endif else begin
            print,'Running FITEXY for ',gals[i]
            fitexy_result = mpfitexy(xdat, ydat, xerr, yerr, guess=[3d,0d],$
               /fixint, /reduce, scatter=fitexy_scat, e_int_reduce=e_int_reduce, $
               /quiet, /silent)
            ymod = xdat*fitexy_result[0]+fitexy_result[1]
            npts = n_elements(xdat)
            boot_ind = boot_indices(npts, nsample = nsamp)
            marr = dblarr(nsamp)
            barr = dblarr(nsamp)
            scatarr = dblarr(nsamp)
            for j=0,nsamp-1 do begin
               boot_result = mpfitexy(xdat[boot_ind[j,*]],ydat[boot_ind[j,*]],$
                  xerr[boot_ind[j,*]],yerr[boot_ind[j,*]],guess=[3d,0d],$
                  /fixint, /reduce, scatter=boot_scat, e_int_reduce=e_int_reduce, $
                  /quiet, /silent)
               marr[j] = boot_result[0]
               barr[j] = boot_result[1]
               scatarr[j] = boot_scat
            endfor
            fitexy_merr = stddev(marr-fitexy_result[0])
            fitexy_berr = stddev(barr-fitexy_result[1])
            fitexy_scaterr = stddev(scatarr-fitexy_scat)
            fitexypar.result[i,*] = fitexy_result
            fitexypar.nfit[i] = n_elements(ymod)
            fitexypar.ymod[i,0:n_elements(ymod)-1] = ymod
            fitexypar.scat[i] = fitexy_scat
            fitexypar.merr[i] = fitexy_merr
            fitexypar.berr[i] = fitexy_berr
            fitexypar.scaterr[i] = fitexy_scaterr
         endelse
         ; Flag if fit is significant; use r_xy not consistent with 0 at 95.5% conf.
         issig = 1b
         if corr[0] eq 0 OR $
            (corr[0] lt 0 AND corr[0]+corr[4] ge 0) OR $
            (corr[0] gt 0 AND corr[0]-corr[3] le 0) then $
            issig = 0b
         cgplot,xran,[0,0],/xsty,/ysty,xran=xran,yran=yran,aspect=1d,$
            pos=pos[*,i],noerase=noerase,title=plotstr.name
         if issig then $
            cgpolygon,[xdat,reverse(xdat),xdat[0]],$
                   [yfit2lo,reverse(yfit2hi),yfit2lo[0]],$
                   color='Medium Gray',fcolor='Medium Gray',thick=1,/fill,noclip=0
         cgoplot,x,y,psym=3,/err_clip,err_color='Slate Gray',err_width=0d,$
                 err_xlow=xerrlo,err_xhi=xerrhi,$
                 err_ylow=yerrlo,err_yhigh=yerrhi
         cgoplot,x,y,psym=16,symsize=0.5
         ; cgoplot,x,yboth,psym=16,symsize=0.5,color='Cyan'
         if issig then begin
            cgoplot,xdat,yfit,color='BLUE',thick=8
            cgoplot,xdat,yfit+sqrt(sigsqr[0]),color='BLUE',thick=4,/linesty
            cgoplot,xdat,yfit-sqrt(sigsqr[0]),color='BLUE',thick=4,/linesty
            ; fitexy model
            cgoplot,xdat,ymod,color='Magenta',thick=4
            cgoplot,xdat,ymod+fitexy_scat,color='Magenta',thick=4,/linesty
            cgoplot,xdat,ymod-fitexy_scat,color='Magenta',thick=4,/linesty
         endif
         if linmixpar.pval[i,1] eq 1b then pvalll='<' else pvalll=''
         printf,lun_tmp,plotstr.name,amp,linmixpar.nfit[i],amp,$
;            corr[0],'$_{-',corr[1],'}^{+',corr[2],'}$ ',amp,$
            corr[0],' (',pvalll,linmixpar.pval[i,0],')',amp,$
            alpha[0],'$_{-',alpha[1],'}^{+',alpha[2],'}$ ',amp,$
            beta[0],'$_{-',beta[1],'}^{+',beta[2],'}$ ',amp,$
            sqrt(sigsqr[0]),'$_{-',sigsqr[1]/(2d*sqrt(sigsqr[0])),'}^{+',$
            sigsqr[2]/(2d*sqrt(sigsqr[0])),'}$ ',dslash,$
            format='(A10,A3,I4,A3,'+$
;            'D6.2,A0,D6.2,A0,D6.2,A0,A3,'+$
            'D6.2,A0,A0,G0.1,A0,A3,'+$
            'D6.2,A0,D6.2,A0,D6.2,A0,A3,'+$
            'D6.2,A0,D6.2,A0,D6.2,A0,A3,'+$
            'D6.2,A0,D6.2,A0,D6.2,A0,A3)'
         printf,lun_tmp,'',amp,'',amp,'',amp,$
            '0.0',amp,$
            fitexy_result[0],'$\pm$',fitexy_merr,amp,$
            fitexy_scat,'$\pm$',fitexy_scaterr,dslash,$
            format='(A10,A3,A4,A3,A6,A3,A0,A3,'+$
            'D6.2,A0,D0.2,A3,D6.2,A0,D0.2,A3)'
         fit1_rxyall[i]=corr[0]
         fit1_weq0all[i]=alpha[0]
         fit1_mall_linmix[i]=beta[0]
         fit1_mall_fitexy[i]=fitexy_result[0]
         fit1_dweqall_linmix[i]=sqrt(sigsqr[0])
         fit1_dweqall_fitexy[i]=fitexy_scat
;;        from Richmond et al. 1994
;         ebv_modr = dindgen(11)/10d*0.4d
;         weq_modr = -0.0956d + 4.91d*ebv_modr
;         cgoplot,ebv_modr,weq_modr,color='Cyan'
;;        relationship from Poznanski et al. 2012; scatter is negligible (80 mA)
;         ebv_modp = dindgen(101)/100d*2d
;         logebv_modp = alog10(ebv_modp)
;         weq_modp = 1/1.17d*(logebv_modp+1.85d)
;         cgoplot,ebv_modp,weq_modp,color='Black',thick=4
;;        relationship from Baron et al. 2016
;         ebv_modb = dindgen(11)/10d*0.3d
;         weq_modb = 10.08d*ebv_modb
;         cgoplot,ebv_modb,weq_modb,color=105,thick=4
      endif
   endfor
   xtit='Stellar E(B-V)'
   ytit='Total '+textoidl('W_{eq}')+' ($\angstrom$)'
   cgtext,xtit,0.5d,0.02d,align=0.5,/norm,chars=1.5
   cgtext,ytit,0.03d,0.5d,align=0.5,orient=90d,/norm,chars=1.5
   
   if redolinmix_totweq_v_ebv_stel then $
      save,linmixpar,file=plotdir+'totweq_v_ebv_stel_linmix.xdr'
   if redofitexy_totweq_v_ebv_stel then $
      save,fitexypar,file=plotdir+'totweq_v_ebv_stel_fitexy.xdr'

   cgps_close

   free_lun,lun_tmp
   

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; NaD: Weq vs. Stellar E(B-V)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
     
   nx = 4
   ny = 2
   cgps_open,plotdir+'s7nad_weq_v_ebv_stel.eps',charsize=1,$
             /encap,/nomatch,/inches,xs=plotquantum*nx,ys=plotquantum*ny,/qui
;  units are pixels per tenths of an inch
   pos = cglayout([nx,ny],oxmar=[4,0],oymar=[4,3],ixmar=[0,0],$
                  unit=!D.X_PX_CM*2.54d/10d,xgap=4d,ygap=5d)
   noerase=0b
   for i=0,7 do begin
      xdr = mapdir+gals[i]+'/vormerge/'+gals[i]+'_plots.xdr'
      if file_test(xdr) then begin
         restore,file=xdr
         x = plotstr.stel_ebv
         xerrlo = plotstr.stel_errebv[*,*,0]
         xerrhi = plotstr.stel_errebv[*,*,1]
         y1 = plotstr.weq_abs_and_em
         yerrlo1 = plotstr.weqerrlo_abs_and_em
         yerrhi1 = plotstr.weqerrhi_abs_and_em
         y2 = plotstr.weq_em_and_abs
         yerrlo2 = plotstr.weqerrlo_em_and_abs
         yerrhi2 = plotstr.weqerrhi_em_and_abs
         y3 = plotstr.weq_abs_or_em
         yerrlo3 = plotstr.weqerrlo_abs_or_em
         yerrhi3 = plotstr.weqerrhi_abs_or_em
         y4 = plotstr.weq_noabs_noem
         yerrlo4 = plotstr.weqerrlo_noabs_noem
         yerrhi4 = plotstr.weqerrhi_noabs_noem
         yall = [y1,y2,y3,y4]
         xall = [x,x,x,x]
         igd = where(xall ne bad AND yall ne bad)
         xran = [0.95d,1.05d]*[min(xall[igd]),max(xall[igd])]
         if xran[0] le 0d then xran[0] = -0.05d*max(xall[igd])
         yran = [min(yall[igd]),max(yall[igd])]
         if yran[0] lt 0 then yran *= [1.05d,1.05d] $
         else if yran[0] eq 0 then yran = [-0.05d*yran[1],1.05*yran[1]] $
         else yran *= [0.95d,1.05d]
         if i gt 0 then noerase=1b
         cgplot,xran,[0,0],/xsty,/ysty,xran=xran,yran=yran,aspect=1d,$
                pos=pos[*,i],noerase=noerase
         cgoplot,x,y1,psym=3,/err_clip,err_color='Gray',err_width=0d,$
                err_xlow=xerrlo,err_xhi=xerrhi,$
                err_ylow=yerrlo1,err_yhigh=yerrhi1
         cgoplot,x,y2,psym=3,/err_clip,err_color='Gray',err_width=0d,$
                err_xlow=xerrlo,err_xhi=xerrhi,$
                err_ylow=yerrlo2,err_yhigh=yerrhi2
         cgoplot,x,y3,psym=3,/err_clip,err_color='Gray',err_width=0d,$
                err_xlow=xerrlo,err_xhi=xerrhi,$
                err_ylow=yerrlo3,err_yhigh=yerrhi3
         cgoplot,x,y4,psym=3,/err_clip,err_color='Gray',err_width=0d,$
                err_xlow=xerrlo,err_xhi=xerrhi,$
                err_ylow=yerrlo4,err_yhigh=yerrhi4
         cgoplot,x,y1,psym=16,symsize=0.5,color='Blue'
         cgoplot,x,y2,psym=16,symsize=0.5,color='Red'
         cgoplot,x,y3,psym=16,symsize=0.5
         cgoplot,x,y4,psym=16,symsize=0.5
      endif
   endfor
   xtit='stellar E(B-V)'
   ytit=textoidl('W_{eq}')+' ($\angstrom$)'
   
   cgps_close


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; NaD: Total Weq vs. Gas E(B-V)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
     
   nx = 4
   ny = 2
   ngal = nx*ny
   
   if keyword_set(redolinmix) then redolinmix_totweq_v_ebv_gas = 1b $
   else redolinmix_totweq_v_ebv_gas = 0b
   ; redolinmix_totweq_v_ebv_gas = 1b
   if keyword_set(redofitexy) then redofitexy_totweq_v_ebv_gas = 1b $
   else redofitexy_totweq_v_ebv_gas = 0b
   ; redofitexy_totweq_v_ebv_gas = 1b
    
   openw,lun_tmp,tabdir+'table_totweq_v_ebv_gas.txt',/get_lun

   cgps_open,plotdir+'s7nad_totweq_v_ebv_gas.eps',charsize=1,$
             /encap,/nomatch,/inches,xs=plotquantum*nx*0.95d,ys=plotquantum*ny,/qui
;  units are pixels per tenths of an inch
   pos = cglayout([nx,ny],oxmar=[6,0],oymar=[6,3],ixmar=[-1,0],$
                  unit=!D.X_PX_CM*2.54d/10d,xgap=4d,ygap=5d)
   noerase=0b
   fit2_rxyall=dblarr(ngal)
   fit2_weq0all=dblarr(ngal)
   fit2_mall_linmix=dblarr(ngal)
   fit2_mall_fitexy=dblarr(ngal)
   fit2_dweqall_linmix=dblarr(ngal)
   fit2_dweqall_fitexy=dblarr(ngal)
   maxpts = 1000
   if redolinmix_totweq_v_ebv_gas then begin
      linmixpar = {alpha:dblarr(ngal,5),$
                   beta:dblarr(ngal,5),$
                   corr:dblarr(ngal,5),$
                   sigsqr:dblarr(ngal,5),$
                   nfit:dblarr(ngal),$
                   pval:dblarr(ngal,2),$
                   xdat:dblarr(ngal,maxpts)+bad,$
                   ydat:dblarr(ngal,maxpts)+bad,$
                   xerr:dblarr(ngal,maxpts)+bad,$
                   yerr:dblarr(ngal,maxpts)+bad,$
                   yfit:dblarr(ngal,maxpts)+bad,$
                   yfitlo:dblarr(ngal,maxpts)+bad,$
                   yfithi:dblarr(ngal,maxpts)+bad,$
                   yfit2lo:dblarr(ngal,maxpts)+bad,$
                   yfit2hi:dblarr(ngal,maxpts)+bad}
   endif else begin
      fitfile = plotdir+'totweq_v_ebv_gas_linmix.xdr'
      if file_test(fitfile) then $
         restore,file=fitfile $
      else stop,'WARNING: File '+fitfile+' does not exist; aborting.'
   endelse
   if redofitexy_totweq_v_ebv_gas then begin
      fitexypar = {result:dblarr(ngal,2),$
         nfit:intarr(ngal),$
         ymod:dblarr(ngal,maxpts)+bad,$
         scat:dblarr(ngal),$
         berr:dblarr(ngal),$
         merr:dblarr(ngal),$
         scaterr:dblarr(ngal)}
   endif else begin
      fitfile = plotdir+'totweq_v_ebv_gas_fitexy.xdr'
      if file_test(fitfile) then $
         restore,file=fitfile $
      else stop,'WARNING: File '+fitfile+' does not exist; aborting.'
   endelse
   for i=0,ngal-1 do begin
      xdr = mapdir+gals[i]+'/vormerge/'+gals[i]+'_plots.xdr'
      if file_test(xdr) then begin
         restore,file=xdr
         x = plotstr.ebv
         xerrlo = plotstr.errebv[*,*]
         xerrhi = plotstr.errebv[*,*]
         y = plotstr.weq_abs_plus_em
         yerrlo = plotstr.weqerrlo_abs_plus_em
         yerrhi = plotstr.weqerrhi_abs_plus_em
         igd = where(x ne bad AND y ne bad)
         xran = [-0.049,2.199]
         yran = [-4,11]
         if i gt 0 then noerase=1b
         if ~ redolinmix_totweq_v_ebv_gas then begin
            sigsqr = reform(linmixpar.sigsqr[i,*],5)
            alpha = reform(linmixpar.alpha[i,*],5)
            beta = reform(linmixpar.beta[i,*],5)
            corr = reform(linmixpar.corr[i,*],5)
            xdat = reform(linmixpar.xdat[i,0:linmixpar.nfit[i]-1],$
                          linmixpar.nfit[i])
            ydat = reform(linmixpar.ydat[i,0:linmixpar.nfit[i]-1],$
                          linmixpar.nfit[i])
            xerr = reform(linmixpar.xerr[i,0:linmixpar.nfit[i]-1],$
                          linmixpar.nfit[i])
            yerr = reform(linmixpar.yerr[i,0:linmixpar.nfit[i]-1],$
                          linmixpar.nfit[i])
            yfit = reform(linmixpar.yfit[i,0:linmixpar.nfit[i]-1],$
                          linmixpar.nfit[i])
            yfitlo = reform(linmixpar.yfitlo[i,0:linmixpar.nfit[i]-1],$
                            linmixpar.nfit[i])
            yfithi = reform(linmixpar.yfithi[i,0:linmixpar.nfit[i]-1],$
                            linmixpar.nfit[i])
            yfit2lo = reform(linmixpar.yfit2lo[i,0:linmixpar.nfit[i]-1],$
                            linmixpar.nfit[i])
            yfit2hi = reform(linmixpar.yfit2hi[i,0:linmixpar.nfit[i]-1],$
                            linmixpar.nfit[i])
         endif else begin
            print,'Running LINMIX for ',gals[i]
            ifit = where(x ne bad AND y ne bad)
;           get rid of duplicates due to Voronoi tiling
            xdatall = x[ifit]
            xerrall = (xerrlo[ifit]+xerrhi[ifit])/2d
            ydatall = y[ifit]
            yerrall = (yerrlo[ifit]+yerrhi[ifit])/2d
            zdatall = xdatall+ydatall
;           Unique elements from x+y
            isort = bsort(zdatall)           
            iuniq = uniq(zdatall[isort])
            isortuniq = isort[iuniq]
            xdat_unsort = xdatall[isortuniq]
            ydat_unsort = ydatall[isortuniq]
            xerr_unsort = xerrall[isortuniq]
            yerr_unsort = yerrall[isortuniq]
            iresort = bsort(xdat_unsort)
            xdat = xdat_unsort[iresort]
            ydat = ydat_unsort[iresort]
            xerr = xerr_unsort[iresort]
            yerr = yerr_unsort[iresort]
            fitout = drt_runlinmix(xdat,ydat,xerr=xerr,yerr=yerr,/metro,$
                                   alpha=alpha,beta=beta,corr=corr,$
                                   sigsqr=sigsqr,pval=pval,miniter=lm_miniter)
            print,pval
            yfit = fitout[*,0]
            yfitlo = fitout[*,1]
            yfithi = fitout[*,2]
            yfit2lo = fitout[*,3]
            yfit2hi = fitout[*,4]
            nfit = n_elements(xdat)
            linmixpar.alpha[i,*] = alpha
            linmixpar.beta[i,*] = beta
            linmixpar.corr[i,*] = corr
            linmixpar.sigsqr[i,*] = sigsqr
            linmixpar.nfit[i] = nfit
            linmixpar.pval[i,*] = pval
            linmixpar.xdat[i,0:nfit-1] = xdat
            linmixpar.ydat[i,0:nfit-1] = ydat
            linmixpar.xerr[i,0:nfit-1] = xerr
            linmixpar.yerr[i,0:nfit-1] = yerr
            linmixpar.yfit[i,0:nfit-1] = yfit
            linmixpar.yfitlo[i,0:nfit-1] = yfitlo
            linmixpar.yfithi[i,0:nfit-1] = yfithi
            linmixpar.yfit2lo[i,0:nfit-1] = yfit2lo
            linmixpar.yfit2hi[i,0:nfit-1] = yfit2hi
         endelse
         ; FITEXY fits
         if ~ redofitexy_totweq_v_ebv_stel then begin
            fitexy_result = reform(fitexypar.result[i,*],2)
            nfit = fitexypar.nfit[i]
            ymod = reform(fitexypar.ymod[i,0:nfit-1],nfit)
            fitexy_scat = fitexypar.scat[i]
            fitexy_merr = fitexypar.merr[i]
            fitexy_berr = fitexypar.berr[i]
            fitexy_scaterr = fitexypar.scaterr[i]
         endif else begin
            print,'Running FITEXY for ',gals[i]
            fitexy_result = mpfitexy(xdat, ydat, xerr, yerr, guess=[3d,0d],$
               /fixint, /reduce, scatter=fitexy_scat, e_int_reduce=e_int_reduce, $
               /quiet, /silent)
            ymod = xdat*fitexy_result[0]+fitexy_result[1]
            npts = n_elements(xdat)
            boot_ind = boot_indices(npts, nsample = nsamp)
            marr = dblarr(nsamp)
            barr = dblarr(nsamp)
            scatarr = dblarr(nsamp)
            for j=0,nsamp-1 do begin
               boot_result = mpfitexy(xdat[boot_ind[j,*]],ydat[boot_ind[j,*]],$
                  xerr[boot_ind[j,*]],yerr[boot_ind[j,*]],guess=[3d,0d],$
                  /fixint, /reduce, scatter=boot_scat, e_int_reduce=e_int_reduce, $
                  /quiet, /silent)
               marr[j] = boot_result[0]
               barr[j] = boot_result[1]
               scatarr[j] = boot_scat
            endfor
            fitexy_merr = stddev(marr-fitexy_result[0])
            fitexy_berr = stddev(barr-fitexy_result[1])
            fitexy_scaterr = stddev(scatarr-fitexy_scat)
            fitexypar.result[i,*] = fitexy_result
            fitexypar.nfit[i] = n_elements(ymod)
            fitexypar.ymod[i,0:n_elements(ymod)-1] = ymod
            fitexypar.scat[i] = fitexy_scat
            fitexypar.merr[i] = fitexy_merr
            fitexypar.berr[i] = fitexy_berr
            fitexypar.scaterr[i] = fitexy_scaterr
         endelse
         ; Flag if fit is significant; use r_xy not consistent with 0 at 95.5% conf.
         issig = 1b
         if corr[0] eq 0 OR $
            (corr[0] lt 0 AND corr[0]+corr[4] ge 0) OR $
            (corr[0] gt 0 AND corr[0]-corr[3] le 0) then $
            issig = 0b
         cgplot,xran,[0,0],/xsty,/ysty,xran=xran,yran=yran,aspect=1d,$
                pos=pos[*,i],noerase=noerase,title=plotstr.name
         if issig then $
            cgpolygon,[xdat,reverse(xdat),xdat[0]],$
                   [yfit2lo,reverse(yfit2hi),yfit2lo[0]],$
                   color='Medium Gray',fcolor='Medium Gray',thick=1,/fill,noclip=0
         cgoplot,x,y,psym=3,/err_clip,err_color='Slate Gray',err_width=0d,$
                 err_xlow=xerrlo,err_xhi=xerrhi,$
                 err_ylow=yerrlo,err_yhigh=yerrhi
         cgoplot,x,y,psym=16,symsize=0.5
         if issig then begin
            cgoplot,xdat,yfit,color='BLUE',thick=8
            cgoplot,xdat,yfit+sqrt(sigsqr[0]),color='BLUE',thick=4,/linesty
            cgoplot,xdat,yfit-sqrt(sigsqr[0]),color='BLUE',thick=4,/linesty
            ; fitexy model
            cgoplot,xdat,ymod,color='Magenta',thick=4
            cgoplot,xdat,ymod+fitexy_scat,color='Magenta',/linesty,thick=4
            cgoplot,xdat,ymod-fitexy_scat,color='Magenta',/linesty,thick=4
         endif
         if linmixpar.pval[i,1] eq 1b then pvalll='<' else pvalll=''
         printf,lun_tmp,plotstr.name,amp,linmixpar.nfit[i],amp,$
;            corr[0],'$_{-',corr[1],'}^{+',corr[2],'}$ ',amp,$
            corr[0],' (',pvalll,linmixpar.pval[i,0],')',amp,$
            alpha[0],'$_{-',alpha[1],'}^{+',alpha[2],'}$ ',amp,$
            beta[0],'$_{-',beta[1],'}^{+',beta[2],'}$ ',amp,$
            sqrt(sigsqr[0]),'$_{-',sigsqr[1]/(2d*sqrt(sigsqr[0])),'}^{+',$
            sigsqr[2]/(2d*sqrt(sigsqr[0])),'}$ ',dslash,$
            format='(A10,A3,I4,A3,'+$
 ;           'D6.2,A0,D6.2,A0,D6.2,A0,A3,'+$
            'D6.2,A0,A0,G0.1,A0,A3,'+$
            'D6.2,A0,D6.2,A0,D6.2,A0,A3,'+$
            'D6.2,A0,D6.2,A0,D6.2,A0,A3,'+$
            'D6.2,A0,D6.2,A0,D6.2,A0,A3)'
         printf,lun_tmp,'',amp,'',amp,$
            '0.0',amp,$
            fitexy_result[0],'$\pm$',fitexy_merr,amp,$
            fitexy_scat,'$\pm$',fitexy_scaterr,dslash,$
            format='(A10,A3,A4,A3,A6,A3,'+$
            'D6.2,A0,D0.2,A3,D6.2,A0,D0.2,A3)'
         fit2_rxyall[i]=corr[0]
         fit2_weq0all[i]=alpha[0]
         fit2_mall_linmix[i]=beta[0]
         fit2_mall_fitexy[i]=fitexy_result[0]
         fit2_dweqall_linmix[i]=sqrt(sigsqr[0])
         fit2_dweqall_fitexy[i]=fitexy_scat
;;        from Richmond et al. 1994
;         ebv_modr = dindgen(11)/10d*0.4d
;         weq_modr = -0.0956d + 4.91d*ebv_modr
;         cgoplot,ebv_modr,weq_modr,color='Cyan'
;;        relationship from Poznanski et al. 2012; scatter is negligible (80 mA)
;         ebv_modp = dindgen(101)/100d*2d
;         logebv_modp = alog10(ebv_modp)
;         weq_modp = 1/1.17d*(logebv_modp+1.85d)
;         cgoplot,ebv_modp,weq_modp,color='Cyan'
;;        relationship from Baron et al. 2016
;         ebv_modb = dindgen(11)/10d*0.3d
;         weq_modb = 10.08d*ebv_modb
;         cgoplot,ebv_modb,weq_modb,color='Cyan'


      endif
   endfor
   xtit='Gas E(B-V)'
   ytit='Total '+textoidl('W_{eq}')+' ($\angstrom$)'
   cgtext,xtit,0.5d,0.02d,align=0.5,/norm,chars=1.5
   cgtext,ytit,0.03d,0.5d,align=0.5,orient=90d,/norm,chars=1.5
   
   if redolinmix_totweq_v_ebv_gas then $
      save,linmixpar,file=plotdir+'totweq_v_ebv_gas_linmix.xdr'
   if redofitexy_totweq_v_ebv_gas then $
      save,fitexypar,file=plotdir+'totweq_v_ebv_gas_fitexy.xdr'

   cgps_close

   free_lun,lun_tmp
   

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; NaD: Weq vs. Gas E(B-V)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;     
;   nx = 4
;   ny = 2
;   cgps_open,plotdir+'s7nad_weq_v_ebv_gas.eps',charsize=1,$
;             /encap,/nomatch,/inches,xs=plotquantum*nx,ys=plotquantum*ny,/qui
;;  units are pixels per tenths of an inch
;   pos = cglayout([nx,ny],oxmar=[4,0],oymar=[4,3],ixmar=[0,0],$
;                  unit=!D.X_PX_CM*2.54d/10d,xgap=4d,ygap=5d)
;   noerase=0b
;   for i=0,7 do begin
;      xdr = mapdir+gals[i]+'/vormerge/'+gals[i]+'_plots.xdr'
;      if file_test(xdr) then begin
;         restore,file=xdr
;         x = plotstr.ebv
;         xerrlo = plotstr.errebv[*,*]
;         xerrhi = plotstr.errebv[*,*]
;         y1 = plotstr.weq_abs_and_em
;         yerrlo1 = plotstr.weqerrlo_abs_and_em
;         yerrhi1 = plotstr.weqerrhi_abs_and_em
;         y2 = plotstr.weq_em_and_abs
;         yerrlo2 = plotstr.weqerrlo_em_and_abs
;         yerrhi2 = plotstr.weqerrhi_em_and_abs
;         y3 = plotstr.weq_abs_or_em
;         yerrlo3 = plotstr.weqerrlo_abs_or_em
;         yerrhi3 = plotstr.weqerrhi_abs_or_em
;         y4 = plotstr.weq_noabs_noem
;         yerrlo4 = plotstr.weqerrlo_noabs_noem
;         yerrhi4 = plotstr.weqerrhi_noabs_noem
;         yall = [y1,y2,y3,y4]
;         xall = [x,x,x,x]
;         igd = where(xall ne bad AND yall ne bad)
;         xran = [0.95d,1.05d]*[min(xall[igd]),max(xall[igd])]
;         if xran[0] le 0d then xran[0] = -0.05d*max(xall[igd])
;         yran = [min(yall[igd]),max(yall[igd])]
;         if yran[0] lt 0 then yran *= [1.05d,1.05d] $
;         else if yran[0] eq 0 then yran = [-0.05d*yran[1],1.05*yran[1]] $
;         else yran *= [0.95d,1.05d]
;         if i gt 0 then noerase=1b
;         cgplot,xran,[0,0],/xsty,/ysty,xran=xran,yran=yran,aspect=1d,$
;                pos=pos[*,i],noerase=noerase
;         cgoplot,x,y1,psym=3,/err_clip,err_color='Gray',err_width=0d,$
;                err_xlow=xerrlo,err_xhi=xerrhi,$
;                err_ylow=yerrlo1,err_yhigh=yerrhi1
;         cgoplot,x,y2,psym=3,/err_clip,err_color='Gray',err_width=0d,$
;                err_xlow=xerrlo,err_xhi=xerrhi,$
;                err_ylow=yerrlo2,err_yhigh=yerrhi2
;         cgoplot,x,y3,psym=3,/err_clip,err_color='Gray',err_width=0d,$
;                err_xlow=xerrlo,err_xhi=xerrhi,$
;                err_ylow=yerrlo3,err_yhigh=yerrhi3
;         cgoplot,x,y4,psym=3,/err_clip,err_color='Gray',err_width=0d,$
;                err_xlow=xerrlo,err_xhi=xerrhi,$
;                err_ylow=yerrlo4,err_yhigh=yerrhi4
;         cgoplot,x,y1,psym=16,symsize=0.5,color='Blue'
;         cgoplot,x,y2,psym=16,symsize=0.5,color='Red'
;         cgoplot,x,y3,psym=16,symsize=0.5
;         cgoplot,x,y4,psym=16,symsize=0.5
;      endif
;   endfor
;   xtit='Gas E(B-V)'
;   ytit=textoidl('W_{eq}')+' ($\angstrom$)'
;   
;   cgps_close

   
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; NaD: Weq vs. Gas E(B-V)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
     
   nx = 4
   ny = 2
   ngal = nx*ny
    
   xran = [-1.3d,0.3d]
   yran = [-1d,1d]

;   if keyword_set(redolinmix) then redolinmix_weq_v_ebv_gas = 1b $
;   else redolinmix_weq_v_ebv_gas = 0b

   weq_v_ebv_gas_x = []
   weq_v_ebv_gas_y = []
   weq_v_ebv_gas_perarea_x = []
   weq_v_ebv_gas_perarea_y = []
   weq_v_ebv_gas_perbin_x = []
   weq_v_ebv_gas_perbin_y = []

;   openw,lun_weq_v_ebv,tabdir+'table_weq_v_ebv_gas.txt',/get_lun

   cgps_open,plotdir+'s7nad_weq_v_ebv_gas.eps',charsize=1,$
             /encap,/nomatch,/inches,xs=plotquantum*nx*0.95d,ys=plotquantum*ny,/qui
;  units are pixels per tenths of an inch
   pos = cglayout([nx,ny],oxmar=[6,0],oymar=[6,3],ixmar=[-1,0],$
                  unit=!D.X_PX_CM*2.54d/10d,xgap=4d,ygap=5d)
   noerase=0b
;   if redolinmix_weq_v_ebv_gas then begin
;      maxpts = 1000
;      linmixpar = {alpha:dblarr(ngal,3),$
;                   beta:dblarr(ngal,3),$
;                   corr:dblarr(ngal,3),$
;                   sigsqr:dblarr(ngal,3),$
;                   nfit:dblarr(ngal),$
;                   yfit:dblarr(ngal,maxpts)+bad,$
;                   yfitlo:dblarr(ngal,maxpts)+bad,$
;                   yfithi:dblarr(ngal,maxpts)+bad}
;   endif else begin
;      fitfile = plotdir+'weq_v_ebv_gas_fits.xdr'
;      if file_test(fitfile) then $
;         restore,file=fitfile $
;      else stop,'WARNING: File '+fitfile+' does not exist; aborting.'
;   endelse
   for i=0,ngal-1 do begin
      xdr = mapdir+gals[i]+'/vormerge/'+gals[i]+'_plots.xdr'
      if file_test(xdr) then begin
         restore,file=xdr
         x = plotstr.ebv
         xerrlo = plotstr.errebv[*,*]
         xerrhi = plotstr.errebv[*,*]
         y = plotstr.weq_abs
         yerrlo = plotstr.weqerrlo_abs
         yerrhi = plotstr.weqerrhi_abs
         igd = where(x ne bad AND x ne 0 AND y ne bad AND y ne 0)
         xlin = x
         xlinerrlo = xerrlo
         xlinerrhi = xerrhi
         ylin = y
         ylinerrlo = yerrlo
         ylinerrhi = yerrhi
         x = alog10(x[igd])
         y = alog10(y[igd])
         weq_v_ebv_gas_x = [weq_v_ebv_gas_x,x]
         weq_v_ebv_gas_y = [weq_v_ebv_gas_y,y]
         for j=0,intareanorm[i]-1 do begin
            weq_v_ebv_gas_perarea_x = [weq_v_ebv_gas_perarea_x,x]
            weq_v_ebv_gas_perarea_y = [weq_v_ebv_gas_perarea_y,y]
         endfor
         xerrlo = x - alog10(xlin[igd] - xerrlo[igd])
         xerrhi = alog10(xlin[igd] + xerrhi[igd]) - x
         yerrlo = y - alog10(ylin[igd] - yerrlo[igd])
         yerrhi = alog10(ylin[igd] + yerrhi[igd]) - y
;        fix NaNs due to very large low errors
         inanx = where(~ finite(xerrlo),ctnanx)
         if ctnanx gt 0 then xerrlo[inanx] = xerrhi[inanx]
         inany = where(~ finite(yerrlo),ctnany)
         if ctnany gt 0 then yerrlo[inany] = yerrhi[inany]
;        get rid of duplicates due to Voronoi tiling
         ifit = where(x ne 0d AND x ne bad AND y ne bad)
         xfitall = x[ifit]
         xerrfitall = (xerrlo[ifit]+xerrhi[ifit])/2d
         yfitall = y[ifit]
         yerrfitall = (yerrlo[ifit]+yerrhi[ifit])/2d
         zfitall = xfitall+yfitall
;        Unique elements from x+y
         isort = bsort(zfitall)
         iuniq = uniq(zfitall[isort])
         isortuniq = isort[iuniq]
         xfit_unsort = xfitall[isortuniq]
         yfit_unsort = yfitall[isortuniq]
         xfiterr_unsort = xerrfitall[isortuniq]
         yfiterr_unsort = yerrfitall[isortuniq]
         iresort = bsort(xfit_unsort)
         xfit = xfit_unsort[iresort]
         yfit = yfit_unsort[iresort]
         xfiterr = xfiterr_unsort[iresort]
         yfiterr = yfiterr_unsort[iresort]
         ;
         weq_v_ebv_gas_perbin_x = [weq_v_ebv_gas_perbin_x,xfit]
         weq_v_ebv_gas_perbin_y = [weq_v_ebv_gas_perbin_y,yfit]
         if i gt 0 then noerase=1b
;         if ~ redolinmix_weq_v_ebv_gas then begin
;            sigsqr = reform(linmixpar.sigsqr[i,*],3)
;            alpha = reform(linmixpar.alpha[i,*],3)
;            beta = reform(linmixpar.beta[i,*],3)
;            corr = reform(linmixpar.corr[i,*],3)
;            yfitout = reform(linmixpar.yfit[i,0:linmixpar.nfit[i]-1],$
;                             linmixpar.nfit[i])
;            yfitoutlo = reform(linmixpar.yfitlo[i,0:linmixpar.nfit[i]-1],$
;                             linmixpar.nfit[i])
;            yfitouthi = reform(linmixpar.yfithi[i,0:linmixpar.nfit[i]-1],$
;                             linmixpar.nfit[i])
;         endif else begin
;            print,'Running LINMIX for ',gals[i]
;            fitout = drt_runlinmix(xfit,yfit,xerr=xfiterr,yerr=yfiterr,$
;                                   alpha=alpha,beta=beta,corr=corr,sigsqr=sigsqr)
;            yfitout = fitout[*,0]
;            yfitoutlo = fitout[*,1]
;            yfitouthi = fitout[*,2]
;            nfit = n_elements(xfit)
;            linmixpar.alpha[i,*] = alpha
;            linmixpar.beta[i,*] = beta
;            linmixpar.corr[i,*] = corr
;            linmixpar.sigsqr[i,*] = sigsqr
;            linmixpar.nfit[i] = nfit
;            linmixpar.yfit[i,0:nfit-1] = yfitout
;            linmixpar.yfitlo[i,0:nfit-1] = yfitoutlo
;            linmixpar.yfithi[i,0:nfit-1] = yfitouthi
;         endelse
         cgplot,[0],/xsty,/ysty,xran=xran,yran=yran,aspect=1d,$
                pos=pos[*,i],noerase=noerase,title=plotstr.name,/nodata
;         cgpolygon,[xfit,reverse(xfit),xfit[0]],$
;                   [yfitoutlo,reverse(yfitouthi),yfitoutlo[0]],$
;                   color='grey',fcolor='grey',thick=1,/fill,noclip=0
         cgoplot,x,y,psym=3,/err_clip,err_color='Slate Gray',err_width=0d,$
                 err_xlow=xerrlo,err_xhi=xerrhi,$
                 err_ylow=yerrlo,err_yhigh=yerrhi
         cgoplot,x,y,psym=16,symsize=0.5
;         cgoplot,xfit,yfitout,color='BLUE',thick=8
;         cgoplot,xfit,yfitout+sqrt(sigsqr[1]),color='BLUE',thick=4,/linesty
;         cgoplot,xfit,yfitout-sqrt(sigsqr[1]),color='BLUE',thick=4,/linesty
;         printf,lun_weq_v_ebv,gals[i],amp,$
;                linmixpar.nfit[i],amp,$
;                corr[1],'$_{-',corr[1]-corr[0],'}^{+',corr[2]-corr[1],'}$ ',amp,$
;                alpha[1],'$_{-',alpha[1]-alpha[0],'}^{+',alpha[2]-alpha[1],'}$ ',amp,$
;                beta[1],'$_{-',beta[1]-beta[0],'}^{+',beta[2]-beta[1],'}$ ',amp,$
;                sqrt(sigsqr[1]),'$_{-',sqrt(sigsqr[1])-sqrt(sigsqr[0]),'}^{+',sqrt(sigsqr[2])-sqrt(sigsqr[1]),'}$ ',dslash,$
;                format='(A10,A3,I4,A3,'+$
;                       'D6.2,A0,D6.2,A0,D6.2,A0,A3,'+$
;                       'D6.2,A0,D6.2,A0,D6.2,A0,A3,'+$
;                       'D6.2,A0,D6.2,A0,D6.2,A0,A3,'+$
;                       'D6.2,A0,D6.2,A0,D6.2,A0,A3)'
;;        from Richmond et al. 1994
;         ebv_modr = dindgen(11)/10d*0.4d
;         weq_modr = -0.0956d + 4.91d*ebv_modr
;         cgoplot,alog10(ebv_modr),alog10(weq_modr),color='Cyan'
;        relationship from Poznanski et al. 2012; scatter is negligible (80 mA)
         ebv_modp = dindgen(101)/100d*2d
         logebv_modp = alog10(ebv_modp)
         weq_modp = 1/1.17d*(logebv_modp+1.85d)
         cgoplot,logebv_modp,alog10(weq_modp),color='Black',thick=4
         cgoplot,logebv_modp,alog10(weq_modp+0.080d),color='Black',thick=2,linesty=2
         cgoplot,logebv_modp,alog10(weq_modp-0.080d),color='Black',thick=2,linesty=2
;        relationship from Baron et al. 2016
         ebv_modb = dindgen(11)/10d*0.3d
         weq_modb = 10.08d*ebv_modb
         cgoplot,alog10(ebv_modb),alog10(weq_modb),color=105,thick=4
         cgoplot,alog10(ebv_modb),alog10(weq_modb+0.55d),color=105,thick=2,linesty=2
         cgoplot,alog10(ebv_modb),alog10(weq_modb-0.55d),color=105,thick=2,linesty=2
;;        relationship from Chen et al. 2010, estimated
;         ebv_modb = dindgen(11)/10d*(2.5d - 1.3d)/3.1d + 1.3d/3.1d
;         weq_modb = 1d/(1.3d/3.1d)*ebv_modb
;         cgoplot,alog10(ebv_modb),alog10(weq_modb),color=102,thick=8
;
;
;         ;  LIRGs
;         cgoplot,logebv_lirgs,alog10(weq_lirgs),psym=15,color=103,symsize=1
;         cgoplot,logebv_lirgs,alog10(weq_lirgs),psym=6,symsize=1


; data range
cgoplot,[-1.7,0.3],[0.1,1.05],color=101,thick=8,linesty=1
cgoplot,[-1.7,0.3],[-1.4,0.3],color=101,thick=8,linesty=1


      endif
   endfor
   xtit='log [ Gas E(B-V) ]'
   ytit='log [ W$\downeq$ / $\angstrom$ ]'
   cgtext,xtit,0.5d,0.02d,align=0.5,/norm,chars=1.5
   cgtext,ytit,0.03d,0.5d,align=0.5,orient=90d,/norm,chars=1.5
   
;   if redolinmix_weq_v_ebv_gas then $
;      save,linmixpar,file=plotdir+'weq_v_ebv_gas_fits.xdr'

   cgps_close

  binsizex = 0.2d
  binsizey = 0.2d
  density = hist_2d(weq_v_ebv_gas_x,weq_v_ebv_gas_y,$
                    bin1=binsizex,bin2=binsizey,$
                    min1=xran[0],min2=yran[0],max1=xran[1],max2=yran[1])
  dsize = size(density)
  density = double(density[0:dsize[1]-2,0:dsize[2]-2])
  maxDensity = double(max(density)) ; Ceil(Max(density)/1d2) * 1d2
  normDensity = density / maxDensity
  scaledDensity = BytScl(density, Min=0, Max=0.5*maxDensity,top=150)
  inz = where(scaledDensity ne 0)
  scaledDensity[inz]+=42
  density = hist_2d(weq_v_ebv_gas_perarea_x,weq_v_ebv_gas_perarea_y,$
                    bin1=binsizex,bin2=binsizey,$
                    min1=xran[0],min2=yran[0],max1=xran[1],max2=yran[1])
  density = double(density[0:dsize[1]-2,0:dsize[2]-2])
  maxDensity = double(max(density)) ;Ceil(Max(density)/1d2) * 1d2
  normDensity_perarea = density / maxDensity
  scaledDensity_perarea = BytScl(density, Min=0, Max=maxDensity,top=150)
  inz = where(scaledDensity_perarea ne 0)
  scaledDensity_perarea[inz]+=10
  density = hist_2d(weq_v_ebv_gas_perbin_x,weq_v_ebv_gas_perbin_y,$
                    bin1=binsizex,bin2=binsizey,$
                    min1=xran[0],min2=yran[0],max1=xran[1],max2=yran[1])
  density = double(density[0:dsize[1]-2,0:dsize[2]-2])
  maxDensity = double(max(density)) ;Ceil(Max(density)/1d2) * 1d2
  normDensity_perbin = density / maxDensity
  scaledDensity_perbin = BytScl(density, Min=0, Max=0.75*maxDensity,top=150)
  inz = where(scaledDensity_perbin ne 0)
  scaledDensity_perbin[inz]+=10

  cgps_open,plotdir+'s7nad_weq_v_ebv_gas_all.eps',charsize=1,$
            /encap,/nomatch,/inches,xs=plotquantum*3.5d,ys=plotquantum*2d,/qui
  xcon = (dindgen(dsize[1]-1)/double(dsize[1]-2))*(xran[1]-xran[0]-binsizex)+xran[0]+binsizex/2d
  ycon = (dindgen(dsize[2]-1)/double(dsize[2]-2))*(yran[1]-yran[0]-binsizey)+yran[0]+binsizey/2d
  for i=0,1 do begin
     cgloadct,0,/reverse
     if i eq 0 then begin
        cgimage,rebin(scaledDensity_perbin,(dsize[1]-1)*samplefac2,$
                      (dsize[2]-1)*samplefac2,/sample),$
                xran=xran,yran=yran,xtit=xtit,ytit=ytit,$
                Position=[0.075, 0.125, 0.475, 0.95], /Axes, /keep
        tvlct,[[27],[158],[119]],100
        tvlct,[[217],[95],[2]],101
        tvlct,[[117],[112],[179]],102
        tvlct,[[231],[41],[138]],103
        tvlct,[[0],[0],[0]],104
        tvlct,[[102],[166],[30]],105
        ;  cgplot,[0],/xsty,/ysty,xran=xran,yran=yran,aspect=1d,/nodata,xtit=xtit,ytit=ytit
        ;        relationship from Poznanski et al. 2012; scatter is negligible (80 mA)
        ebv_modp = dindgen(101)/100d*2d
        logebv_modp = alog10(ebv_modp)
        weq_modp = 1/1.17d*(logebv_modp+1.85d)
        cgoplot,logebv_modp,alog10(weq_modp),color='Black',thick=4
        cgoplot,logebv_modp,alog10(weq_modp+0.080d),color='Black',thick=2,linesty=2
        cgoplot,logebv_modp,alog10(weq_modp-0.080d),color='Black',thick=2,linesty=2
        ;          relationship from Baron et al. 2016
        ebv_modb = dindgen(11)/10d*0.3d
        weq_modb = 10.08d*ebv_modb
        cgoplot,alog10(ebv_modb),alog10(weq_modb),color=105,thick=4
        cgoplot,alog10(ebv_modb),alog10(weq_modb+0.55d),color=105,thick=2,linesty=2
        cgoplot,alog10(ebv_modb),alog10(weq_modb-0.55d),color=105,thick=2,linesty=2
        ;  relationship from Chen et al. 2010, estimated
        ebv_modb = dindgen(11)/10d*(2.5d - 1.3d)/3.1d + 1.3d/3.1d
        weq_modb = 1d/(1.3d/3.1d)*ebv_modb
        cgoplot,alog10(ebv_modb),alog10(weq_modb),color=102,thick=8        
        ; data range
        cgoplot,[-1.7,0.3],[0.1,1.05],color=101,thick=8,linesty=1
        cgoplot,[-1.7,0.3],[-1.4,0.3],color=101,thick=8,linesty=1


        cgoplot,weq_v_ebv_gas_perbin_x,weq_v_ebv_gas_perbin_y,psym=16,symsize=0.2
        cgcontour,normDensity_perbin,xcon,ycon,/noerase,/over,c_color='Red',$
                  levels=[0.2,0.4,0.8],c_thick=4
        maxbin = max(normDensity_perbin,pixmax1d)
        pixmax = array_indices(normDensity_perbin,pixmax1d)
        weq_v_ebv_gas_perbin_max = dblarr(2)
        weq_v_ebv_gas_perbin_max[0] = $
          interpol(xcon,dindgen(dsize[1]-1),[pixmax[0]])
        weq_v_ebv_gas_perbin_max[1] = $
          interpol(ycon,dindgen(dsize[2]-1),[pixmax[1]])
        cgoplot,[weq_v_ebv_gas_perbin_max[0]],$
                [weq_v_ebv_gas_perbin_max[1]],psym=9,symsize=2,color='BLU3',$
                err_xlow=binsizex/2d,err_xhi=binsizex/2d,$
                err_ylow=binsizey/2d,err_yhi=binsizey/2d,$
                err_thick=4d,thick=4d
        pixcent = centroid(normDensity_perbin)
        weq_v_ebv_gas_perbin_centroid = dblarr(2)
        weq_v_ebv_gas_perbin_centroid[0] = $
          interpol(xcon,dindgen(dsize[1]-1),[pixcent[0]])
        weq_v_ebv_gas_perbin_centroid[1] = $
          interpol(ycon,dindgen(dsize[2]-1),[pixcent[1]])
        cgoplot,[weq_v_ebv_gas_perbin_centroid[0]],$
                [weq_v_ebv_gas_perbin_centroid[1]],$
                psym=16,symsize=2,color='RED3',$
                err_xlow=binsizex/2d,err_xhi=binsizex/2d,$
                err_ylow=binsizey/2d,err_yhi=binsizey/2d,$
                err_thick=4d
        pklab = string('Peak=[',weq_v_ebv_gas_perbin_max[0],',',$
                       weq_v_ebv_gas_perbin_max[1],']',$
                       format='(A0,D0.1,A0,D0.1,A0)')
        comlab = string('CofM=[',weq_v_ebv_gas_perbin_centroid[0],',',$
                        weq_v_ebv_gas_perbin_centroid[1],']',$
                        format='(A0,D0.2,A0,D0.2,A0)')
        al_legend,[pklab,comlab],col=['BLU3','RED3'],psym=[9,16],symsize=[2,2],$
                  /left,/top,box=0,spacing=1.5,thick=[4d,4d]
        cgtext,'Unweighted',-0.7,-0.9,align=0.5
;        result = gauss2dfit(scaledDensity_perbin,a,xcon,ycon,/tilt)
;        cgoplot,[a[4]],[a[5]],psym=16,symsize=1,color='Orange'
;        cgoplot,[a[4]+a[2]*cos(a[6]),a[4]-a[2]*cos(a[6])],$
;                [a[5]-a[2]*sin(a[6]),a[5]+a[2]*sin(a[6])],$
;                color='Orange',thick=4
;        cgoplot,[a[4]-a[3]*sin(a[6]),a[4]+a[3]*sin(a[6])],$
;                [a[5]-a[3]*cos(a[6]),a[5]+a[3]*cos(a[6])],$
;                color='Orange',thick=4
     endif else begin
        cgimage,rebin(scaledDensity_perarea,(dsize[1]-1)*samplefac2,$
                      (dsize[2]-1)*samplefac2,/sample),xran=xran,yran=yran,$
                Position=[0.55, 0.125, 0.95, 0.95], /Axes, /noerase, /keep

        tvlct,[[27],[158],[119]],100
        tvlct,[[217],[95],[2]],101
        tvlct,[[117],[112],[179]],102
        tvlct,[[231],[41],[138]],103
        tvlct,[[0],[0],[0]],104
        tvlct,[[102],[166],[30]],105
        ;  cgplot,[0],/xsty,/ysty,xran=xran,yran=yran,aspect=1d,/nodata,xtit=xtit,ytit=ytit
        ;        relationship from Poznanski et al. 2012; scatter is negligible (80 mA)
        ebv_modp = dindgen(101)/100d*2d
        logebv_modp = alog10(ebv_modp)
        weq_modp = 1/1.17d*(logebv_modp+1.85d)
        cgoplot,logebv_modp,alog10(weq_modp),color='Black',thick=4
        cgoplot,logebv_modp,alog10(weq_modp+0.080d),color='Black',thick=2,linesty=2
        cgoplot,logebv_modp,alog10(weq_modp-0.080d),color='Black',thick=2,linesty=2
        ;          relationship from Baron et al. 2016
        ebv_modb = dindgen(11)/10d*0.3d
        weq_modb = 10.08d*ebv_modb
        cgoplot,alog10(ebv_modb),alog10(weq_modb),color=105,thick=4
        cgoplot,alog10(ebv_modb),alog10(weq_modb+0.55d),color=105,thick=2,linesty=2
        cgoplot,alog10(ebv_modb),alog10(weq_modb-0.55d),color=105,thick=2,linesty=2
        ;  relationship from Chen et al. 2010, estimated
        ebv_modb = dindgen(11)/10d*(2.5d - 1.3d)/3.1d + 1.3d/3.1d
        weq_modb = 1d/(1.3d/3.1d)*ebv_modb
        cgoplot,alog10(ebv_modb),alog10(weq_modb),color=102,thick=8        
        ; data range
        cgoplot,[-1.7,0.3],[0.1,1.05],color=101,thick=8,linesty=1
        cgoplot,[-1.7,0.3],[-1.4,0.3],color=101,thick=8,linesty=1

        cgoplot,weq_v_ebv_gas_perbin_x,weq_v_ebv_gas_perbin_y,psym=16,symsize=0.2
        cgcontour,normDensity_perarea,xcon,ycon,/noerase,/over,c_color='Red',$
                  levels=[0.2,0.4,0.8],c_thick=4
        maxbin = max(normDensity_perarea,pixmax1d)
        pixmax = array_indices(normDensity_perarea,pixmax1d)
        weq_v_ebv_gas_perarea_max = dblarr(2)
        weq_v_ebv_gas_perarea_max[0] = $
           interpol(xcon,dindgen(dsize[1]-1),[pixmax[0]])
        weq_v_ebv_gas_perarea_max[1] = $
           interpol(ycon,dindgen(dsize[2]-1),[pixmax[1]])
        cgoplot,[weq_v_ebv_gas_perarea_max[0]],$
                [weq_v_ebv_gas_perarea_max[1]],$
                psym=9,symsize=2,color='BLU3',$
                err_xlow=binsizex/2d,err_xhi=binsizex/2d,$
                err_ylow=binsizey/2d,err_yhi=binsizey/2d,$
                err_thick=4d,thick=4d
        pixcent = centroid(normDensity_perarea)
        weq_v_ebv_gas_perarea_centroid = dblarr(2)
        weq_v_ebv_gas_perarea_centroid[0] = $
          interpol(xcon,dindgen(dsize[1]-1),[pixcent[0]])
        weq_v_ebv_gas_perarea_centroid[1] = $
          interpol(ycon,dindgen(dsize[2]-1),[pixcent[1]])
        cgoplot,[weq_v_ebv_gas_perarea_centroid[0]],$
                [weq_v_ebv_gas_perarea_centroid[1]],$
                psym=16,symsize=2,color='RED3',$
                err_xlow=binsizex/2d,err_xhi=binsizex/2d,$
                err_ylow=binsizey/2d,err_yhi=binsizey/2d,$
                err_thick=4d
        pklab = string('Peak=[',weq_v_ebv_gas_perarea_max[0],',',$
                       weq_v_ebv_gas_perarea_max[1],']',$
                       format='(A0,D0.1,A0,D0.1,A0)')
        comlab = string('CofM=[',weq_v_ebv_gas_perarea_centroid[0],',',$
                        weq_v_ebv_gas_perarea_centroid[1],']',$
                        format='(A0,D0.2,A0,D0.2,A0)')
        al_legend,[pklab,comlab],col=['BLU3','RED3'],psym=[9,16],symsize=[2,2],$
                  /left,/top,box=0,spacing=1.5,thick=[4d,4d]
        cgtext,'Area-weighted',-0.7,-0.9,align=0.5
;        result = gauss2dfit(scaledDensity_perarea,a,xcon,ycon,/tilt)
;        cgoplot,[a[4]],[a[5]],psym=16,symsize=1,color='Orange'
;        cgoplot,[a[4]+a[2]*cos(a[6]),a[4]-a[2]*cos(a[6])],$
;                [a[5]-a[2]*sin(a[6]),a[5]+a[2]*sin(a[6])],$
;                color='Orange',thick=4
;        cgoplot,[a[4]-a[3]*sin(a[6]),a[4]+a[3]*sin(a[6])],$
;                [a[5]-a[3]*cos(a[6]),a[5]+a[3]*cos(a[6])],$
;                color='Orange',thick=4
     endelse
  
  
  
endfor

  cgps_close

   
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; NaD: Total Weq vs. Stellar color; g-i colors only
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
     
   nx = 4
   ny = 2

   if keyword_set(redolinmix) then redolinmix_totweq_v_col = 1b $
   else redolinmix_totweq_v_col = 0b
   ; redolinmix_totweq_v_col = 1b
   if keyword_set(redofitexy) then redofitexy_totweq_v_col = 1b $
   else redofitexy_totweq_v_col = 0b
   ; redofitexy_totweq_v_col = 1b

   openw,lun_tmp,tabdir+'table_totweq_v_col.txt',/get_lun

   cgps_open,plotdir+'s7nad_totweq_v_col.eps',charsize=1,$
             /encap,/nomatch,/inches,xs=plotquantum*nx*0.95d,ys=plotquantum*ny,/qui
;  units are pixels per tenths of an inch
   pos = cglayout([nx,ny],oxmar=[6,0],oymar=[6,3],ixmar=[-1,0],$
                  unit=!D.X_PX_CM*2.54d/10d,xgap=4d,ygap=5d)
   fit3_rxyall=dblarr(ngal)
   fit3_giall_linmix=dblarr(ngal)
   fit3_giall_fitexy=dblarr(ngal)
   fit3_mall_linmix=dblarr(ngal)
   fit3_mall_fitexy=dblarr(ngal)
   fit3_dweqall_linmix=dblarr(ngal)
   fit3_dweqall_fitexy=dblarr(ngal)
   maxpts = 1000
   if redolinmix_totweq_v_col then begin
      linmixpar = {alpha:dblarr(ngal,5),$
                   beta:dblarr(ngal,5),$
                   corr:dblarr(ngal,5),$
                   sigsqr:dblarr(ngal,5),$
                   nfit:dblarr(ngal),$
                   pval:dblarr(ngal,2),$
                   xdat:dblarr(ngal,maxpts)+bad,$
                   ydat:dblarr(ngal,maxpts)+bad,$
                   xerr:dblarr(ngal,maxpts)+bad,$
                   yerr:dblarr(ngal,maxpts)+bad,$
                   yfit:dblarr(ngal,maxpts)+bad,$
                   yfitlo:dblarr(ngal,maxpts)+bad,$
                   yfithi:dblarr(ngal,maxpts)+bad,$
                   yfit2lo:dblarr(ngal,maxpts)+bad,$
                   yfit2hi:dblarr(ngal,maxpts)+bad}
   endif else begin
      fitfile = plotdir+'totweq_v_col_linmix.xdr'
      if file_test(fitfile) then $
         restore,file=fitfile $
      else stop,'WARNING: File '+fitfile+' does not exist; aborting.'
   endelse
   if redofitexy_totweq_v_col then begin
      fitexypar = {result:dblarr(ngal,2),$
         nfit:intarr(ngal),$
         ymod:dblarr(ngal,maxpts)+bad,$
         scat:dblarr(ngal),$
         berr:dblarr(ngal),$
         merr:dblarr(ngal),$
         scaterr:dblarr(ngal)}
   endif else begin
      fitfile = plotdir+'totweq_v_col_fitexy.xdr'
      if file_test(fitfile) then $
         restore,file=fitfile $
      else stop,'WARNING: File '+fitfile+' does not exist; aborting.'
   endelse
   noerase=0b
   i=0
   j=0
   gals_gicol = ['ngc1266','ngc1808','eso500','ngc5728','eso339','ic5063','ic5169','ic1481']
   for k=0,ngal-1 do begin
      xdr = mapdir+gals[k]+'/vormerge/'+gals[k]+'_plots.xdr'
      if file_test(xdr) AND where(gals[k] eq gals_gicol) ne -1 then begin
         restore,file=xdr
         x = plotstr.cmap
         if n_elements(x) gt 1 then begin
            xerrlo = plotstr.ecmap
            xerrhi = plotstr.ecmap
            y = plotstr.weq_abs_plus_em
            yboth = plotstr.weq_abs_and_em
            yerrlo = plotstr.weqerrlo_abs_plus_em
            yerrhi = plotstr.weqerrhi_abs_plus_em
            igd = where(x ne bad AND y ne bad)
;            xran = [0.95d,1.05d]*[min(x[igd]),max(x[igd])]
;            if xran[0] le 0d then xran[0] = -0.05d*max(x[igd])
;            yran = [min(y[igd]),max(y[igd])]
;            if yran[0] lt 0 then yran *= [1.05d,1.05d] $
;            else if yran[0] eq 0 then yran = [-0.05d*yran[1],1.05*yran[1]] $
;            else yran *= [0.95d,1.05d]
;            xran = [xran[0],xran[0]+1.2]
            xran = [0.4,1.8]
            if k eq 1 then xran = [1.3,3.3]
            yran = [-4,10]
            if i gt 0 then noerase=1b
            if ~ redolinmix_totweq_v_col then begin
               sigsqr = reform(linmixpar.sigsqr[i,*],5)
               alpha = reform(linmixpar.alpha[i,*],5)
               beta = reform(linmixpar.beta[i,*],5)
               corr = reform(linmixpar.corr[i,*],5)
               xdat = reform(linmixpar.xdat[i,0:linmixpar.nfit[i]-1],$
                             linmixpar.nfit[i])
               ydat = reform(linmixpar.ydat[i,0:linmixpar.nfit[i]-1],$
                             linmixpar.nfit[i])
               xerr = reform(linmixpar.xerr[i,0:linmixpar.nfit[i]-1],$
                             linmixpar.nfit[i])
               yerr = reform(linmixpar.yerr[i,0:linmixpar.nfit[i]-1],$
                             linmixpar.nfit[i])
               yfit = reform(linmixpar.yfit[i,0:linmixpar.nfit[i]-1],$
                             linmixpar.nfit[i])
               yfitlo = reform(linmixpar.yfitlo[i,0:linmixpar.nfit[i]-1],$
                               linmixpar.nfit[i])
               yfithi = reform(linmixpar.yfithi[i,0:linmixpar.nfit[i]-1],$
                                linmixpar.nfit[i])
               yfit2lo = reform(linmixpar.yfit2lo[i,0:linmixpar.nfit[i]-1],$
                                linmixpar.nfit[i])
               yfit2hi = reform(linmixpar.yfit2hi[i,0:linmixpar.nfit[i]-1],$
                                linmixpar.nfit[i])
            endif else begin
               print,'Running LINMIX for ',gals[i]
               ifit = where(x ne 0d AND x ne bad AND y ne bad)
;              get rid of duplicates due to Voronoi tiling
               xdatall = x[ifit]
               xerrall = (xerrlo[ifit]+xerrhi[ifit])/2d
               ydatall = y[ifit]
               yerrall = (yerrlo[ifit]+yerrhi[ifit])/2d
               zdatall = xdatall+ydatall
;              Unique elements from x+y
               isort = bsort(zdatall)           
               iuniq = uniq(zdatall[isort])
               isortuniq = isort[iuniq]
               xdat_unsort = xdatall[isortuniq]
               ydat_unsort = ydatall[isortuniq]
               xerr_unsort = xerrall[isortuniq]
               yerr_unsort = yerrall[isortuniq]
               iresort = bsort(xdat_unsort)
               xdat = xdat_unsort[iresort]
               ydat = ydat_unsort[iresort]
               xerr = xerr_unsort[iresort]
               yerr = yerr_unsort[iresort]
               fitout = drt_runlinmix(xdat,ydat,xerr=xerr,yerr=yerr,/metro,$
                                   alpha=alpha,beta=beta,corr=corr,$
                                   sigsqr=sigsqr,pval=pval,miniter=lm_miniter)
               print,pval
               yfit = fitout[*,0]
               yfitlo = fitout[*,1]
               yfithi = fitout[*,2]
               yfit2lo = fitout[*,3]
               yfit2hi = fitout[*,4]
               nfit = n_elements(xdat)
               linmixpar.alpha[i,*] = alpha
               linmixpar.beta[i,*] = beta
               linmixpar.corr[i,*] = corr
               linmixpar.sigsqr[i,*] = sigsqr
               linmixpar.nfit[i] = nfit
               linmixpar.pval[i,*] = pval
               linmixpar.xdat[i,0:nfit-1] = xdat
               linmixpar.ydat[i,0:nfit-1] = ydat
               linmixpar.xerr[i,0:nfit-1] = xerr
               linmixpar.yerr[i,0:nfit-1] = yerr
               linmixpar.yfit[i,0:nfit-1] = yfit
               linmixpar.yfitlo[i,0:nfit-1] = yfitlo
               linmixpar.yfithi[i,0:nfit-1] = yfithi
               linmixpar.yfit2lo[i,0:nfit-1] = yfit2lo
               linmixpar.yfit2hi[i,0:nfit-1] = yfit2hi
            endelse
            ; FITEXY fits
            if ~ redofitexy_totweq_v_col then begin
               fitexy_result = reform(fitexypar.result[i,*],2)
               nfit = fitexypar.nfit[i]
               ymod = reform(fitexypar.ymod[i,0:nfit-1],nfit)
               fitexy_scat = fitexypar.scat[i]
               fitexy_merr = fitexypar.merr[i]
               fitexy_berr = fitexypar.berr[i]
               fitexy_scaterr = fitexypar.scaterr[i]
            endif else begin
               print,'Running FITEXY for ',gals[i]
               fitexy_result = mpfitexy(xdat, ydat, xerr, yerr, guess=[3d,0d],$
                  /reduce, scatter=fitexy_rms, e_int_reduce=fitexy_scat, $
                  /quiet, /silent)
               ymod = xdat*fitexy_result[0]+fitexy_result[1]
               npts = n_elements(xdat)
               boot_ind = boot_indices(npts, nsample = nsamp)
               marr = dblarr(nsamp)
               barr = dblarr(nsamp)
               scatarr = dblarr(nsamp)
               for l=0,nsamp-1 do begin
                  boot_result = mpfitexy(xdat[boot_ind[l,*]],ydat[boot_ind[l,*]],$
                     xerr[boot_ind[l,*]],yerr[boot_ind[l,*]],guess=[3d,0d],$
                     /reduce, scatter=boot_rms, e_int_reduce=boot_scat, $
                     /quiet, /silent)
                  marr[l] = boot_result[0]
                  barr[l] = boot_result[1]
                  scatarr[l] = boot_scat
               endfor
               fitexy_merr = stddev(marr-fitexy_result[0])
               fitexy_berr = stddev(barr-fitexy_result[1])
               fitexy_scaterr = stddev(scatarr-fitexy_scat)
               fitexypar.result[i,*] = fitexy_result
               fitexypar.nfit[i] = n_elements(ymod)
               fitexypar.ymod[i,0:n_elements(ymod)-1] = ymod
               fitexypar.scat[i] = fitexy_scat
               fitexypar.merr[i] = fitexy_merr
               fitexypar.berr[i] = fitexy_berr
               fitexypar.scaterr[i] = fitexy_scaterr
            endelse
            ; Flag if fit is significant; use r_xy not consistent with 0 at 95.5% conf.
            issig = 1b
            if corr[0] eq 0 OR $
               (corr[0] lt 0 AND corr[0]+corr[4] ge 0) OR $
               (corr[0] gt 0 AND corr[0]-corr[3] le 0) then $
               issig = 0b
            cgplot,xran,[0,0],/xsty,/ysty,xran=xran,yran=yran,aspect=1d,$
                   pos=pos[*,k],noerase=noerase,title=plotstr.name
            if issig then $
               cgpolygon,[xdat,reverse(xdat),xdat[0]],$
                  [yfit2lo,reverse(yfit2hi),yfit2lo[0]],$
                  color='Medium Gray',fcolor='Medium Gray',thick=1,/fill,noclip=0
            cgoplot,x,y,psym=3,/err_clip,err_color='Slate Gray',err_width=0d,$
                    err_ylow=yerrlo,err_yhigh=yerrhi,$
                    err_xlow=xerrlo,err_xhi=xerrhi
            cgoplot,x,y,psym=16,symsize=0.5
            if issig then begin
               cgoplot,xdat,yfit,color='BLUE',thick=8
               cgoplot,xdat,yfit+sqrt(sigsqr[0]),color='BLUE',thick=4,/linesty
               cgoplot,xdat,yfit-sqrt(sigsqr[0]),color='BLUE',thick=4,/linesty
               ; fitexy model
               cgoplot,xdat,ymod,color='Magenta',thick=4
               cgoplot,xdat,ymod+fitexy_scat,color='Magenta',/linesty,thick=4
               cgoplot,xdat,ymod-fitexy_scat,color='Magenta',/linesty,thick=4
            endif
            j++
            xint = alpha[0]/(-beta[0])
            xint_errlo = xint*sqrt((alpha[1]/alpha[0])^2d + $
               (beta[2]/beta[0])^2d)
            xint_errhi = xint*sqrt((alpha[2]/alpha[0])^2d + $
               (beta[1]/beta[0])^2d)
            if linmixpar.pval[i,1] eq 1b then pvalll='<' else pvalll=''
            printf,lun_tmp,plotstr.name,amp,linmixpar.nfit[i],amp,$
;               corr[0],'$_{-',corr[1],'}^{+',corr[2],'}$ ',amp,$
               corr[0],' (',pvalll,linmixpar.pval[i,0],')',amp,$
               xint,'$_{-',xint_errlo,'}^{+',xint_errhi,'}$ ',amp,$
               beta[0],'$_{-',beta[1],'}^{+',beta[2],'}$ ',amp,$
               sqrt(sigsqr[0]),'$_{-',sigsqr[1]/(2d*sqrt(sigsqr[0])),'}^{+',$
               sigsqr[2]/(2d*sqrt(sigsqr[0])),'}$ ',dslash,$
               format='(A10,A3,I4,A3,'+$
;               'D6.2,A0,D6.2,A0,D6.2,A0,A3,'+$
               'D6.2,A0,A0,G0.1,A0,A3,'+$
               'D6.2,A0,D6.2,A0,D6.2,A0,A3,'+$
               'D6.2,A0,D6.2,A0,D6.2,A0,A3,'+$
               'D6.2,A0,D6.2,A0,D6.2,A0,A3)'
            xint = fitexy_result[1]/(-fitexy_result[0])
            xint_err = xint*sqrt((fitexy_berr/fitexy_result[1])^2d + $
               (fitexy_merr/fitexy_result[0])^2d)
            printf,lun_tmp,'',amp,'',amp,'',amp,$
               xint,'$\pm$',xint_err,amp,$
               fitexy_result[0],'$\pm$',fitexy_merr,amp,$
               fitexy_scat,'$\pm$',fitexy_scaterr,dslash,$
               format='(A10,A3,A4,A3,A6,A3,D6.2,A0,D0.2,A3,'+$
               'D6.2,A0,D0.2,A3,D6.2,A0,D0.2,A3)'
            fit3_rxyall[i]=corr[0]
            fit3_giall_linmix[i]=alpha[0]/(-beta[0])
            fit3_giall_fitexy[i]=fitexy_result[1]/(-fitexy_result[0])
            fit3_mall_linmix[i]=beta[0]
            fit3_mall_fitexy[i]=fitexy_result[0]
            fit3_dweqall_linmix[i]=sqrt(sigsqr[0])
            fit3_dweqall_fitexy[i]=fitexy_scat
            i++
         endif
      endif
   endfor
   xtit='Color'
   ytit='Total '+textoidl('W_{eq}')+' ($\angstrom$)'
   cgtext,xtit,0.5d,0.02d,align=0.5,/norm,chars=1.5
   cgtext,ytit,0.04d,0.5d,align=0.5,orient=90d,/norm,chars=1.5

   if redolinmix_totweq_v_col then $
      save,linmixpar,file=plotdir+'totweq_v_col_linmix.xdr'
   if redofitexy_totweq_v_col then $
      save,fitexypar,file=plotdir+'totweq_v_col_fitexy.xdr'
   
   cgps_close

   free_lun,lun_tmp

   
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; NaD: Total Weq vs. Stellar color: Others
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

   nx = 1
   ny = 1

   cgps_open,plotdir+'s7nad_totweq_v_col_others.eps',charsize=1,$
             /encap,/nomatch,/inches,xs=plotquantum*nx*0.9d,ys=plotquantum*ny,/qui
;  units are pixels per tenths of an inch
   pos = cglayout([nx,ny],oxmar=[6,0],oymar=[6,3],ixmar=[-1,0],$
                  unit=!D.X_PX_CM*2.54d/10d,xgap=4d,ygap=5d)
   j=0
   gals_ocol = ['ngc1808']
   ytit=textoidl('W_{eq}')+' ($\angstrom$)'
   for i=0,n_elements(gals_ocol)-1 do begin
      xdr = mapdir+gals_ocol[i]+'/vormerge/'+gals_ocol[i]+'_plots.xdr'
      if file_test(xdr) then begin
         restore,file=xdr
         x = plotstr.cmap
         if n_elements(x) gt 1 then begin
            xerrlo = x*0d
            xerrhi = x*0d
            y = plotstr.weq_abs_plus_em
            yerrlo = plotstr.weqerrlo_abs_plus_em
            yerrhi = plotstr.weqerrhi_abs_plus_em
            igd = where(x ne bad AND y ne bad)
            xran = [0.95d,1.05d]*[min(x[igd]),max(x[igd])]
            if xran[0] le 0d then xran[0] = -0.05d*max(x[igd])
            yran = [min(y[igd]),max(y[igd])]
            if yran[0] lt 0 then yran *= [1.05d,1.05d] $
            else if yran[0] eq 0 then yran = [-0.05d*yran[1],1.05*yran[1]] $
            else yran *= [0.95d,1.05d]
            if i eq 0 then begin
               noerase=0b
               ytit = 'Total '+textoidl('W_{eq}')+' ($\angstrom$)'
               xtit = '(B-I) Color'
            endif else if i eq 1 then begin
               noerase=1b
               ytit = ''
               xtit = '(F547M-F606W) Color'
            endif ;else begin
;               noerase=1b
;               ytit = ''
;               xtit = '(g-i)$\downSkyMapper$ Color'
;            endelse
            cgplot,xran,[0,0],/xsty,/ysty,xran=xran,yran=yran,aspect=1d,$
                   pos=pos[*,j],noerase=noerase,title=plotstr.name,$
                   xtit=xtit,ytit=ytit
            cgoplot,x,y,psym=3,/err_clip,err_color='Gray',err_width=0d,$
                   err_xlow=xerrlo,err_xhi=xerrhi,$
                   err_ylow=yerrlo,err_yhigh=yerrhi
            cgoplot,x,y,psym=16,symsize=0.5
            j++
         endif
      endif
   endfor
   
   cgps_close

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; NaD: Total Weq vs. sigma
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
     
   nx = 4
   ny = 2
   ngal = nx*ny

   weq_tot_sigem_cc = dblarr(ngal) + bad
   weq_tot_sigabs_cc = dblarr(ngal) + bad

   cgps_open,plotdir+'s7nad_totweq_v_sig.eps',charsize=1,$
             /encap,/nomatch,/inches,xs=plotquantum*nx*0.95d,ys=plotquantum*ny,/qui
;  units are pixels per tenths of an inch
   pos = cglayout([nx,ny],oxmar=[6,0],oymar=[6,3],ixmar=[-1,0],$
                  unit=!D.X_PX_CM*2.54d/10d,xgap=4d,ygap=5d)
   noerase=0b
   xran = [-250,250]
   yran = [-4,11]
   for i=0,ngal-1 do begin
      xdr = mapdir+gals[i]+'/vormerge/'+gals[i]+'_plots.xdr'
      if file_test(xdr) then begin
         restore,file=xdr
         xa = plotstr.sig_abs
         xe = plotstr.sig_em
         if n_elements(xa) gt 1 OR n_elements(xe) gt 1 then begin
            y = plotstr.weq_abs_plus_em
            yerrlo = plotstr.weqerrlo_abs_plus_em
            yerrhi = plotstr.weqerrhi_abs_plus_em
            igda = where(xa ne bad and y ne bad)
            igde = where(xe ne bad AND y ne bad)
            xe[igde] *= -1d
            if i gt 0 then noerase=1b
            cgplot,xran,[0,0],/xsty,/ysty,xran=xran,yran=yran,aspect=1d,$
                   pos=pos[*,i],noerase=noerase,title=plotstr.name
            cgoplot,xa,y,psym=3,/err_clip,err_color='Gray',err_width=0d,$
                    err_ylow=yerrlo,err_yhigh=yerrhi
            weq_tot_sigabs_cc[i] = correlate(xa[igda],y[igda])
            if n_elements(xe) gt 1 then begin
               cgoplot,xe,y,psym=3,/err_clip,err_color='Gray',err_width=0d,$
                       err_ylow=yerrlo,err_yhigh=yerrhi
               cgoplot,xe,y,psym=9,symsize=0.5,color='RED3'
               weq_tot_sigem_cc[i] = -1d * correlate(xe[igde],y[igde])
            endif
            cgoplot,xa,y,psym=16,symsize=0.5,color='BLU3'
         endif
      endif
   endfor
   xtit='$\sigma$ (km s$\up-1$)'
   ytit=textoidl('Total W_{eq}')+' ($\angstrom$)'
   cgtext,xtit,0.5d,0.02d,align=0.5,/norm,chars=1.5
   cgtext,ytit,0.03d,0.5d,align=0.5,orient=90d,/norm,chars=1.5

   cgps_close




;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; NaD: N(NaI) vs. color
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
     
   nx = 4
   ny = 2

   cgps_open,plotdir+'s7nad_nnai_v_col.eps',charsize=1,$
             /encap,/nomatch,/inches,xs=plotquantum*nx,ys=plotquantum*ny*1.1d,/qui
;  units are pixels per tenths of an inch
   pos = cglayout([nx,ny],oxmar=[6,0],oymar=[6,3],ixmar=[0,0],$
                  unit=!D.X_PX_CM*2.54d/10d,xgap=4d,ygap=5d)
   noerase=0b
   xran = [0.4,1.799]
   yran = [11.5,15.5]
   gals_nnai = ['ngc1266','eso500','ngc5728','eso339','ic5063','ic5169','ic1481']
   for i=0,6 do begin
      if i gt 0 then noerase=1b
      xdr = mapdir+gals_nnai[i]+'/vormerge/'+gals_nnai[i]+'_plots.xdr'
      restore,file=xdr
      cmap = plotstr.cmap
      ecmap = plotstr.ecmap
      mapabs = plotstr.nnai
      mapabserrlo = plotstr.nnai_el
      mapabserrhi = plotstr.nnai_eh
      igdone = plotstr.igdone_nnai
      istone = plotstr.istone_nnai
      igdtwo = plotstr.igdtwo_nnai
      isttwo = plotstr.isttwo_nnai
;     Compute average points
      xavg = [1.3,1.5,1.7,1.9,2.1,2.3,2.5]
      yavg = dblarr(7) + bad
      ysdev = dblarr(7)
      iall = [igdone,istone]
      if igdtwo[0] ne -1 then begin
         iall = [iall,igdtwo,isttwo]
      endif
      xall = cmap[iall]
      xallerr = ecmap[iall]
      yall = mapabs[iall]
      yallerr = (mapabserrlo[iall]+mapabserrhi[iall])/2d
;     get rid of duplicates due to Voronoi tiling
      zall = xall+yall
      isort = bsort(zall)           
      iuniq = uniq(zall[isort])
      isortuniq = isort[iuniq]
      x_unsort = xall[isortuniq]
      y_unsort = yall[isortuniq]
      iresort = bsort(x_unsort)
      xuniq = x_unsort[iresort]
      yuniq = y_unsort[iresort]
      for j=0,n_elements(xavg)-1 do begin
         ibin = where(xuniq ge xavg[j]-0.1d AND $
                      xuniq lt xavg[j]+0.1d AND $
                      xuniq ne bad,ctbin)
         if ctbin gt 0 then begin
            yavg[j]=mean(yuniq[ibin])
;            yavg[j]=total(yall[ibin]/yallerr[ibin]^2d)/total(1d/yallerr[ibin]^2d)
            ysdev[j]=stddev(yuniq[ibin])
         endif
      endfor
      cgplot,cmap[igdone],mapabs[igdone],/xsty,/ysty,xran=xran,yran=yran,$
             color='Gray',noerase=noerase,pos=pos[*,i],$
             psym=3,aspect=1d,err_width=0,err_color='Gray',/err_clip,$
             err_ylow=mapabserrlo[igdone],err_yhigh=mapabserrhi[igdone],$
             err_xlow=ecmap[igdone],err_xhigh=ecmap[igdone],$
             title=plotstr.name
      cgoplot,cmap[istone],mapabs[istone],color='Gray',psym=3,$
              err_width=0,err_color='Gray',/err_clip,$
              err_ylow=mapabserrlo[istone],err_yhigh=mapabserrhi[istone],$
              err_xlow=ecmap[istone],err_xhigh=ecmap[istone]
      if igdtwo[0] ne -1 then begin
         cgoplot,cmap[igdtwo],mapabs[igdtwo],color='Gray',psym=3,$
                 err_width=0,err_color='Gray',/err_clip,$
                 err_ylow=mapabserrlo[igdtwo],err_yhigh=mapabserrhi[igdtwo],$
                 err_xlow=ecmap[igdtwo],err_xhigh=ecmap[igdtwo]
         cgoplot,cmap[isttwo],mapabs[isttwo],color='Gray',psym=3,$
                 err_width=0,err_color='Gray',/err_clip,$
                 err_ylow=mapabserrlo[isttwo],err_yhigh=mapabserrhi[isttwo],$
                 err_xlow=ecmap[isttwo],err_xhigh=ecmap[isttwo]
      endif
      cgoplot,cmap[igdone],mapabs[igdone],color='Red',psym=9,symsize=0.75d
      if igdtwo[0] ne -1 then $
         cgoplot,cmap[igdtwo],mapabs[igdtwo],psym=16,symsize=0.75d
      plotsym,2,thick=2
      cgoplot,cmap[istone],mapabs[istone],color='Red',psym=8
      if igdtwo[0] ne -1 then $
         cgoplot,cmap[isttwo],mapabs[isttwo],psym=8
      cgoplot,xavg,yavg,psym=6,color='Blue',symsize=2,$
              err_ylow=ysdev,err_yhi=ysdev
   endfor
   xtit='(g-i) Color'
   ytit='Total N(NaI) (cm$\up-2$)'
   cgtext,xtit,0.5d,0.02d,align=0.5,/norm,chars=1.5
   cgtext,ytit,0.03d,0.5d,align=0.5,orient=90d,/norm,chars=1.5
   
   cgps_close

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; NaD: N(NaI) vs. stellar E(B-V)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
     
   nx = 4
   ny = 2

   cgps_open,plotdir+'s7nad_nnai_v_ebv_stel.eps',charsize=1,$
             /encap,/nomatch,/inches,xs=plotquantum*nx,ys=plotquantum*ny*1.1d,/qui
;  units are pixels per tenths of an inch
   pos = cglayout([nx,ny],oxmar=[6,0],oymar=[6,3],ixmar=[-1,0],$
                  unit=!D.X_PX_CM*2.54d/10d,xgap=4d,ygap=5d)
   noerase=0b
;   xran = [-0.05,1.2]
   xran = [-2,0.5]
   yran = [11,15.5]
   for i=0,7 do begin
      if i gt 0 then noerase=1b
      xdr = mapdir+gals[i]+'/vormerge/'+gals[i]+'_plots.xdr'
      restore,file=xdr
      stel_ebv_lin = plotstr.stel_ebv
      stel_ebv = plotstr.stel_ebv
      stel_errebvlo = plotstr.stel_errebv[*,*,0]
      stel_errebvhi = plotstr.stel_errebv[*,*,1]
      mapabs = plotstr.nnai
      mapabserrlo = plotstr.nnai_el
      mapabserrhi = plotstr.nnai_eh
      igdone = plotstr.igdone_nnai
      istone = plotstr.istone_nnai
      igdtwo = plotstr.igdtwo_nnai
      isttwo = plotstr.isttwo_nnai
;     Compute average points
      xavg = [-1.75,-1.25,-0.75,-0.25,0.25]
      yavg = dblarr(5) + bad
      ysdev = dblarr(5)
      iall = [igdone,istone]
      if igdtwo[0] ne -1 then begin
         iall = [iall,igdtwo,isttwo]
      endif
      xall = alog10(stel_ebv[iall])
      yall = mapabs[iall]
      yallerr = (mapabserrlo[iall]+mapabserrhi[iall])/2d
;     get rid of duplicates due to Voronoi tiling
      zall = xall+yall
      isort = bsort(zall)           
      iuniq = uniq(zall[isort])
      isortuniq = isort[iuniq]
      x_unsort = xall[isortuniq]
      y_unsort = yall[isortuniq]
      iresort = bsort(x_unsort)
      xuniq = x_unsort[iresort]
      yuniq = y_unsort[iresort]
      for j=0,n_elements(xavg)-1 do begin
         ibin = where(xuniq ge xavg[j]-0.1d AND $
                      xuniq lt xavg[j]+0.1d AND $
                      xuniq ne bad,ctbin)
         if ctbin gt 1 then begin
            yavg[j]=mean(yuniq[ibin])
;            yavg[j]=total(yall[ibin]/yallerr[ibin]^2d)/total(1d/yallerr[ibin]^2d)
            ysdev[j]=stddev(yuniq[ibin])
         endif
      endfor
      cgplot,alog10(stel_ebv[igdone]),mapabs[igdone],/xsty,/ysty,xran=xran,yran=yran,$
             color='Gray',noerase=noerase,pos=pos[*,i],$
             psym=3,aspect=1d,err_width=0,err_color='Gray',/err_clip,$
             err_ylow=mapabserrlo[igdone],err_yhigh=mapabserrhi[igdone],$
             err_xlow=stel_errebvlo[igdone],err_xhi=stel_errebvhi[igdone],$
             title=plotstr.name
      cgoplot,alog10(stel_ebv[istone]),mapabs[istone],color='Gray',psym=3,$
              err_width=0,err_color='Gray',/err_clip,$
              err_xlow=stel_errebvlo[istone],err_xhi=stel_errebvhi[istone],$
              err_ylow=mapabserrlo[istone],err_yhigh=mapabserrhi[istone]
      if igdtwo[0] ne -1 then begin
         cgoplot,alog10(stel_ebv[igdtwo]),mapabs[igdtwo],color='Gray',psym=3,$
                 err_width=0,err_color='Gray',/err_clip,$
                 err_xlow=stel_errebvlo[igdtwo],err_xhi=stel_errebvhi[igdtwo],$
                 err_ylow=mapabserrlo[igdtwo],err_yhigh=mapabserrhi[igdtwo]
         cgoplot,alog10(stel_ebv[isttwo]),mapabs[isttwo],color='Gray',psym=3,$
                 err_width=0,err_color='Gray',/err_clip,$
                 err_xlow=stel_errebvlo[isttwo],err_xhi=stel_errebvhi[isttwo],$
                 err_ylow=mapabserrlo[isttwo],err_yhigh=mapabserrhi[isttwo]
      endif
      cgoplot,alog10(stel_ebv[igdone]),mapabs[igdone],color='Red',psym=9,symsize=0.75d
      if igdtwo[0] ne -1 then $
         cgoplot,alog10(stel_ebv[igdtwo]),mapabs[igdtwo],psym=16,symsize=0.75d
      plotsym,2,thick=2
      cgoplot,alog10(stel_ebv[istone]),mapabs[istone],color='Red',psym=8
      if igdtwo[0] ne -1 then $
         cgoplot,alog10(stel_ebv[isttwo]),mapabs[isttwo],psym=8
      cgoplot,xavg,yavg,psym=6,color='Blue',symsize=2,$
              err_ylow=ysdev,err_yhi=ysdev

;  Hobbs 1974
;  N(NaI) = 1.7x10^14 E(B-V)^1.8
   ebv_hobbs = dindgen(26)/25d*0.5d
   logebv_hobbs = alog10(ebv_hobbs)
   nnai_hobbs = 14d + alog10(1.7) + 1.8d*logebv_hobbs
   cgoplot,logebv_hobbs,nnai_hobbs,color='Cyan',linesty=2,thick=8
   cgoplot,logebv_hobbs,nnai_hobbs+0.3d,color='Cyan',linesty=2,thick=8
   cgoplot,logebv_hobbs,nnai_hobbs-0.3d,color='Cyan',linesty=2,thick=8

;  Phillips et al. 2013: fit to Sembach+93 + their own MW data
   ebv_sembach_phillips = dindgen(41)/40d*0.4d
   logebv_sembach_phillips = alog10(ebv_sembach_phillips)
   nnai_sembach_phillips = 13.180d + 1.125d*(logebv_sembach_phillips+alog10(3d))
   cgoplot, logebv_sembach_phillips, nnai_sembach_phillips,color='Black',thick=4
   cgoplot, logebv_sembach_phillips, nnai_sembach_phillips+0.26d,color='Black',thick=4
   cgoplot, logebv_sembach_phillips, nnai_sembach_phillips-0.26d,color='Black',thick=4
   cgoplot, logebv_sembach_phillips, nnai_sembach_phillips+1d,color='Black',thick=8,linesty=2

   endfor
   xtit='Stellar E(B-V)'
   ytit='Total N(NaI) (cm$\up-2$)'
   cgtext,xtit,0.5d,0.02d,align=0.5,/norm,chars=1.5
   cgtext,ytit,0.03d,0.5d,align=0.5,orient=90d,/norm,chars=1.5
   
   cgps_close

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; NaD: N(NaI) vs. gas E(B-V)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
     
   nx = 4
   ny = 2

   xran = [-1.3d,0.3d]
   yran = [11d,15.5d]

   nnai_v_ebv_gas_x = []
   nnai_v_ebv_gas_y = []
   nnai_v_ebv_gas_perarea_x = []
   nnai_v_ebv_gas_perarea_y = []
   nnai_v_ebv_gas_perbin_x = []
   nnai_v_ebv_gas_perbin_y = []

   cgps_open,plotdir+'s7nad_nnai_v_ebv_gas.eps',charsize=1,$
             /encap,/nomatch,/inches,xs=plotquantum*nx,ys=plotquantum*ny*1.1d,/qui
;  units are pixels per tenths of an inch
   pos = cglayout([nx,ny],oxmar=[6,0],oymar=[6,3],ixmar=[-1,0],$
                  unit=!D.X_PX_CM*2.54d/10d,xgap=4d,ygap=5d)
   noerase=0b
   for i=0,7 do begin

      xdr = mapdir+gals[i]+'/vormerge/'+gals[i]+'_plots.xdr'
      restore,file=xdr
      
      if i gt 0 then noerase=1b
      cgplot,[0],/nodata,/xsty,/ysty,xran=xran,yran=yran,$
             noerase=noerase,pos=pos[*,i],aspect=1d,title=plotstr.name
;;errors
;;  M31 pionts
;   cgoplot,logebv_m31,nnai_m31,err_xlow=logebv_errlo_m31,err_xhi=logebv_errhi_m31,$
;           err_ylow=nnai_errlo_m31,err_yhi=nnai_errhi_m31,/err_clip,$
;           psym=1,err_color='Grey',err_width=0
;;  MC pionts
;   cgoplot,logebv_lmc,nnai_lmc,err_xlow=logebv_errlo_lmc,err_xhi=logebv_errhi_lmc,$
;           err_ylow=nnai_errlo_lmc,err_yhi=nnai_errhi_lmc,/err_clip,$
;           psym=1,err_color='Grey',err_width=0
;;  SNe points
;   cgoplot,logebv_sne,nnai_sne,err_xlow=logebv_errlo_sne,err_xhi=logebv_errhi_sne,$
;           err_ylow=nnai_err_sne,err_yhi=nnai_err_sne,/err_clip,$
;           psym=1,err_color='Grey',err_width=0
;;data
;;  M31 pionts
;   cgoplot,logebv_m31,nnai_m31,psym=15,color=100,symsize=0.5
;;   cgoplot,logebv_m31,nnai_m31,psym=6,symsize=0.5
;;  MC pionts
;   cgoplot,logebv_lmc,nnai_lmc,psym=15,color=101,symsize=0.5
;;   cgoplot,logebv_lmc,nnai_lmc,psym=6,symsize=0.5
;   cgoplot,alog10(ebv_mc),nnai_mc,psym=15,color=101,symsize=0.5
;;   cgoplot,alog10(ebv_mc),nnai_mc,psym=6,symsize=0.5
;;  SNe points
;   cgoplot,logebv_sne,nnai_sne,psym=15,color=102,symsize=0.5
;;   cgoplot,logebv_sne,nnai_sne,psym=6,symsize=0.5
      
      
      ebv_lin = plotstr.ebv
      ebv = plotstr.ebv
      errebvlo = plotstr.errebv[*,*]
      errebvhi = plotstr.errebv[*,*]
      mapabs = plotstr.nnai
      mapabserrlo = plotstr.nnai_el
      mapabserrhi = plotstr.nnai_eh
      igdone = plotstr.igdone_nnai
      istone = plotstr.istone_nnai
      igdtwo = plotstr.igdtwo_nnai
      isttwo = plotstr.isttwo_nnai
;     Compute average points
      xavg = [-1.25,-1.00,-0.75,-0.50,-0.25,0,0.25,0.50]
      yavg = dblarr(7) + bad
      ysdev = dblarr(7)
      iall = [igdone,istone]
      if igdtwo[0] ne -1 then begin
         iall = [iall,igdtwo,isttwo]
      endif
      xall = alog10(ebv[iall])
      yall = mapabs[iall]
      nnai_v_ebv_gas_x = [nnai_v_ebv_gas_x,xall]
      nnai_v_ebv_gas_y = [nnai_v_ebv_gas_y,yall]
      for j=0,intareanorm[i]-1 do begin
         nnai_v_ebv_gas_perarea_x = [nnai_v_ebv_gas_perarea_x,xall]
         nnai_v_ebv_gas_perarea_y = [nnai_v_ebv_gas_perarea_y,yall]
      endfor
      yallerr = (mapabserrlo[iall]+mapabserrhi[iall])/2d
;     get rid of duplicates due to Voronoi tiling
      zall = xall+yall
      isort = bsort(zall)           
      iuniq = uniq(zall[isort])
      isortuniq = isort[iuniq]
      x_unsort = xall[isortuniq]
      y_unsort = yall[isortuniq]
      iresort = bsort(x_unsort)
      xuniq = x_unsort[iresort]
      yuniq = y_unsort[iresort]
      nnai_v_ebv_gas_perbin_x = [nnai_v_ebv_gas_perbin_x,xuniq]
      nnai_v_ebv_gas_perbin_y = [nnai_v_ebv_gas_perbin_y,yuniq]
      for j=0,n_elements(xavg)-1 do begin
         ibin = where(xuniq ge xavg[j]-0.1d AND $
                      xuniq lt xavg[j]+0.1d AND $
                      xuniq ne bad,ctbin)
         if ctbin gt 1 then begin
            yavg[j]=mean(yuniq[ibin])
;            yavg[j]=total(yall[ibin]/yallerr[ibin]^2d)/total(1d/yallerr[ibin]^2d)
            ysdev[j]=stddev(yuniq[ibin])
         endif
      endfor
      cgoplot,alog10(ebv[igdone]),mapabs[igdone],color='Slate Gray',$
             psym=3,err_width=0,err_color='Slate Gray',/err_clip,$
             err_ylow=mapabserrlo[igdone],err_yhigh=mapabserrhi[igdone],$
             err_xlow=errebvlo[igdone],err_xhi=errebvhi[igdone]
;      izeroebv = where(ebv[igdone] eq 0 AND mapabs[igdone] gt 0,ctzeroebv)
;      if ctzeroebv gt 0 then begin
;         ebvtmp = ebverr[igdone]
;         mapabstmp = mapabs[igdone]
;         mapabserrlotmp = mapabserrlo[igdone]
;         mapabserrhitmp = mapabserrhi[igdone]
;         cgoplot,alog10(ebvtmp[izeroebv]),mapabstmp[izeroebv],$
;                 color='Gray',psym=3,err_width=0,err_color='Gray',/err_clip,$
;                 err_ylow=mapabserrlotmp[izeroebv],$
;                 err_yhigh=mapabserrhitmp[izeroebv]
;         plotsym,6,thick=2
;         cgoplot,alog10(ebvtmp[izeroebv]),mapabstmp[izeroebv],color='Red',psym=8
;      endif
      cgoplot,alog10(ebv[istone]),mapabs[istone],color='Slate Gray',psym=3,$
              err_width=0,err_color='Slate Gray',/err_clip,$
              err_xlow=errebvlo[istone],err_xhi=errebvhi[istone],$
              err_ylow=mapabserrlo[istone],err_yhigh=mapabserrhi[istone]
;      izeroebv = where(ebv[istone] eq 0 AND mapabs[istone] gt 0,ctzeroebv)
;      if ctzeroebv gt 0 then begin
;         print,'yes'
;         ebvtmp = ebverr[istone]
;         mapabstmp = mapabs[istone]
;         mapabserrlotmp = mapabserrlo[istone]
;         mapabserrhitmp = mapabserrhi[istone]
;         cgoplot,alog10(ebvtmp[izeroebv]),mapabstmp[izeroebv],$
;                 color='Gray',psym=3,err_width=0,err_color='Gray',/err_clip,$
;                 err_ylow=mapabserrlotmp[izeroebv],$
;                 err_yhigh=mapabserrhitmp[izeroebv]
;         plotsym,6,thick=2
;         cgoplot,alog10(ebvtmp[izeroebv]),mapabstmp[izeroebv],color='Red',psym=8
;      endif
      if igdtwo[0] ne -1 then begin
         cgoplot,alog10(ebv[igdtwo]),mapabs[igdtwo],color='Slate Gray',psym=3,$
                 err_width=0,err_color='Slate Gray',/err_clip,$
                 err_xlow=errebvlo[igdtwo],err_xhi=errebvhi[igdtwo],$
                 err_ylow=mapabserrlo[igdtwo],err_yhigh=mapabserrhi[igdtwo]
         cgoplot,alog10(ebv[isttwo]),mapabs[isttwo],color='Slate Gray',psym=3,$
                 err_width=0,err_color='Slate Gray',/err_clip,$
                 err_xlow=errebvlo[isttwo],err_xhi=errebvhi[isttwo],$
                 err_ylow=mapabserrlo[isttwo],err_yhigh=mapabserrhi[isttwo]
      endif
      cgoplot,alog10(ebv[igdone]),mapabs[igdone],color='Red',psym=9,symsize=0.75d
      if igdtwo[0] ne -1 then $
         cgoplot,alog10(ebv[igdtwo]),mapabs[igdtwo],psym=16,symsize=0.75d
      plotsym,2,thick=2
      cgoplot,alog10(ebv[istone]),mapabs[istone],color='Red',psym=8
      if igdtwo[0] ne -1 then $
         cgoplot,alog10(ebv[isttwo]),mapabs[isttwo],psym=8
;      cgoplot,xavg,yavg,psym=6,color='Blue',symsize=2,$
;              err_ylow=ysdev,err_yhi=ysdev


;  Hobbs 1974
;  N(NaI) = 1.7x10^14 E(B-V)^1.8
   ebv_hobbs = dindgen(26)/25d*0.5d
   logebv_hobbs = alog10(ebv_hobbs)
   logebv_hobbs_extrap = logebv_hobbs+1d
   nnai_hobbs = 14d + alog10(1.7) + 1.8d*logebv_hobbs
   nnai_hobbs_extrap = 14d + alog10(1.7) + 1.8d*logebv_hobbs_extrap
   cgoplot,logebv_hobbs,nnai_hobbs,color=105,thick=4
   cgoplot,logebv_hobbs,nnai_hobbs+0.3d,color=105,thick=2
   cgoplot,logebv_hobbs,nnai_hobbs-0.3d,color=105,thick=2
   cgoplot,logebv_hobbs_extrap,nnai_hobbs_extrap,color=105,thick=4,linesty=1
   cgoplot,logebv_hobbs_extrap,nnai_hobbs_extrap+0.3d,color=105,thick=2,linesty=1
   cgoplot,logebv_hobbs_extrap,nnai_hobbs_extrap-0.3d,color=105,thick=2,linesty=1

;  Phillips et al. 2013: fit to Sembach+93 + their own MW data
   ebv_sembach_phillips = dindgen(41)/40d*0.4d
   logebv_sembach_phillips = alog10(ebv_sembach_phillips)
   logebv_sembach_phillips_extrap = logebv_sembach_phillips + 1d
   nnai_sembach_phillips = 13.180d + 1.125d*(logebv_sembach_phillips+alog10(3d))
   nnai_sembach_phillips_extrap = 13.180d + 1.125d*(logebv_sembach_phillips_extrap+alog10(3d))
   cgoplot, logebv_sembach_phillips, nnai_sembach_phillips,color='Black',thick=4
   cgoplot, logebv_sembach_phillips, nnai_sembach_phillips+0.26d,color='Black',thick=2
   cgoplot, logebv_sembach_phillips, nnai_sembach_phillips-0.26d,color='Black',thick=2
   cgoplot, logebv_sembach_phillips_extrap, nnai_sembach_phillips_extrap,color='Black',thick=4,linesty=1
   cgoplot, logebv_sembach_phillips_extrap, nnai_sembach_phillips_extrap+0.26d,color='Black',thick=2,linesty=1
   cgoplot, logebv_sembach_phillips_extrap, nnai_sembach_phillips_extrap-0.26d,color='Black',thick=2,linesty=1
;   cgoplot, logebv_sembach_phillips, nnai_sembach_phillips+1d,color='Black',thick=8,linesty=2


cgoplot,[-1.7,0.3],[13.05,15.55],color=101,thick=8,linesty=1
cgoplot,[-1.7,0.3],[10.6,13.35],color=101,thick=8,linesty=1


   endfor
   xtit='log [ Gas E(B-V) ]'
   ytit='log [ N(NaI) / cm$\up-2$ ]'
   cgtext,xtit,0.5d,0.02d,align=0.5,/norm,chars=1.5
   cgtext,ytit,0.03d,0.5d,align=0.5,orient=90d,/norm,chars=1.5
   
   cgps_close


  binsizex = 0.2d
  binsizey = 0.5d
  density = hist_2d(nnai_v_ebv_gas_x,nnai_v_ebv_gas_y,$
                    bin1=binsizex,bin2=binsizey,$
                    min1=xran[0],min2=yran[0],max1=xran[1],max2=yran[1])
  dsize = size(density)
  density = double(density[0:dsize[1]-2,0:dsize[2]-2])
  maxDensity = double(max(density)) ;Ceil(Max(density)/1d2) * 1d2
  normDensity = density / maxDensity
  scaledDensity = BytScl(density, Min=0, Max=0.5*maxDensity,top=150)
  inz = where(scaledDensity ne 0)
  scaledDensity[inz]+=42
  density = hist_2d(nnai_v_ebv_gas_perarea_x,nnai_v_ebv_gas_perarea_y,$
                    bin1=binsizex,bin2=binsizey,$
                    min1=xran[0],min2=yran[0],max1=xran[1],max2=yran[1])
  density = double(density[0:dsize[1]-2,0:dsize[2]-2])
  maxDensity = double(max(density)) ;Ceil(Max(density)/1d2) * 1d2
  normDensity_perarea = density / maxDensity
  scaledDensity_perarea = BytScl(density, Min=0, Max=maxDensity,top=150)
  inz = where(scaledDensity_perarea ne 0)
  scaledDensity_perarea[inz]+=10
  density = hist_2d(nnai_v_ebv_gas_perbin_x,nnai_v_ebv_gas_perbin_y,$
                    bin1=binsizex,bin2=binsizey,$
                    min1=xran[0],min2=yran[0],max1=xran[1],max2=yran[1])
  density = double(density[0:dsize[1]-2,0:dsize[2]-2])
  maxDensity = double(max(density)) ;Ceil(Max(density)/1d2) * 1d2
  normDensity_perbin = density / maxDensity
  scaledDensity_perbin = BytScl(density, Min=0, Max=0.75*maxDensity,top=150)
  inz = where(scaledDensity_perbin ne 0)
  scaledDensity_perbin[inz]+=10
   
   cgps_open,plotdir+'s7nad_nnai_v_ebv_gas_all.eps',charsize=1,$
             /encap,/nomatch,/inches,xs=plotquantum*3.5d,ys=plotquantum*2d,/qui
  xcon = (dindgen(dsize[1]-1)/double(dsize[1]-2))*(xran[1]-xran[0]-binsizex)+xran[0]+binsizex/2d
  ycon = (dindgen(dsize[2]-1)/double(dsize[2]-2))*(yran[1]-yran[0]-binsizey)+yran[0]+binsizey/2d
  for i=0,1 do begin
     cgloadct,0,/reverse
     if i eq 0 then begin
        cgimage,rebin(scaledDensity_perbin,(dsize[1]-1)*samplefac2,$
                      (dsize[2]-1)*samplefac2,/sample),$
                xran=xran,yran=yran,xtit=xtit,ytit=ytit,$
                Position=[0.075, 0.125, 0.475, 0.95], /Axes, /keep

        tvlct,[[27],[158],[119]],100
        tvlct,[[217],[95],[2]],101
        tvlct,[[117],[112],[179]],102
        tvlct,[[231],[41],[138]],103
        tvlct,[[0],[0],[0]],104
        tvlct,[[102],[166],[30]],105
        ;  Hobbs 1974
        ;  N(NaI) = 1.7x10^14 E(B-V)^1.8
        ebv_hobbs = dindgen(26)/25d*0.5d
        logebv_hobbs = alog10(ebv_hobbs)
        logebv_hobbs_extrap = logebv_hobbs+1d
        nnai_hobbs = 14d + alog10(1.7) + 1.8d*logebv_hobbs
        nnai_hobbs_extrap = 14d + alog10(1.7) + 1.8d*logebv_hobbs_extrap
        cgoplot,logebv_hobbs,nnai_hobbs,color=105,thick=4
        cgoplot,logebv_hobbs,nnai_hobbs+0.3d,color=105,thick=2
        cgoplot,logebv_hobbs,nnai_hobbs-0.3d,color=105,thick=2
        cgoplot,logebv_hobbs_extrap,nnai_hobbs_extrap,color=105,thick=4,linesty=1
        cgoplot,logebv_hobbs_extrap,nnai_hobbs_extrap+0.3d,color=105,thick=2,linesty=1
        cgoplot,logebv_hobbs_extrap,nnai_hobbs_extrap-0.3d,color=105,thick=2,linesty=1

        ;  Phillips et al. 2013: fit to Sembach+93 + their own MW data
        ebv_sembach_phillips = dindgen(41)/40d*0.4d
        logebv_sembach_phillips = alog10(ebv_sembach_phillips)
        logebv_sembach_phillips_extrap = logebv_sembach_phillips + 1d
        nnai_sembach_phillips = 13.180d + 1.125d*(logebv_sembach_phillips+alog10(3d))
        nnai_sembach_phillips_extrap = 13.180d + 1.125d*(logebv_sembach_phillips_extrap+alog10(3d))
        cgoplot, logebv_sembach_phillips, nnai_sembach_phillips,color='Black',thick=4
        cgoplot, logebv_sembach_phillips, nnai_sembach_phillips+0.26d,color='Black',thick=2
        cgoplot, logebv_sembach_phillips, nnai_sembach_phillips-0.26d,color='Black',thick=2
        cgoplot, logebv_sembach_phillips_extrap, nnai_sembach_phillips_extrap,color='Black',thick=4,linesty=1
        cgoplot, logebv_sembach_phillips_extrap, nnai_sembach_phillips_extrap+0.26d,color='Black',thick=2,linesty=1
        cgoplot, logebv_sembach_phillips_extrap, nnai_sembach_phillips_extrap-0.26d,color='Black',thick=2,linesty=1
        ;   cgoplot, logebv_sembach_phillips, nnai_sembach_phillips+1d,color='Black',thick=8,linesty=2

        cgoplot,[-1.7,0.3],[13.05,15.55],color=101,thick=8,linesty=1
        cgoplot,[-1.7,0.3],[10.6,13.35],color=101,thick=8,linesty=1

        cgoplot,nnai_v_ebv_gas_perbin_x,nnai_v_ebv_gas_perbin_y,psym=16,symsize=0.2
        cgcontour,normDensity_perbin,xcon,ycon,/noerase,/over,c_color='Red',$
                  levels=[0.2,0.4,0.8],c_thick=4
        maxbin = max(normDensity_perbin,pixmax1d)
        pixmax = array_indices(normDensity_perbin,pixmax1d)
        nnai_v_ebv_gas_perbin_max = dblarr(2)
        nnai_v_ebv_gas_perbin_max[0] = $
          interpol(xcon,dindgen(dsize[1]-1),[pixmax[0]])
        nnai_v_ebv_gas_perbin_max[1] = $
          interpol(ycon,dindgen(dsize[2]-1),[pixmax[1]])
        cgoplot,[nnai_v_ebv_gas_perbin_max[0]],$
                [nnai_v_ebv_gas_perbin_max[1]],psym=9,symsize=2,color='BLU3',$
                err_xlow=binsizex/2d,err_xhi=binsizex/2d,$
                err_ylow=binsizey/2d,err_yhi=binsizey/2d,$
                err_thick=4d,thick=4d
        pixcent = centroid(normDensity_perbin)
        nnai_v_ebv_gas_perbin_centroid = dblarr(2)
        nnai_v_ebv_gas_perbin_centroid[0] = $
          interpol(xcon,dindgen(dsize[1]-1),[pixcent[0]])
        nnai_v_ebv_gas_perbin_centroid[1] = $
          interpol(ycon,dindgen(dsize[2]-1),[pixcent[1]])
        cgoplot,[nnai_v_ebv_gas_perbin_centroid[0]],$
                [nnai_v_ebv_gas_perbin_centroid[1]],$
                psym=16,symsize=2,color='RED3',$
                err_xlow=binsizex/2d,err_xhi=binsizex/2d,$
                err_ylow=binsizey/2d,err_yhi=binsizey/2d,$
                err_thick=4d
        pklab = string('Peak=[',nnai_v_ebv_gas_perbin_max[0],',',$
                       nnai_v_ebv_gas_perbin_max[1],']',$
                       format='(A0,D0.1,A0,D0.1,A0)')
        comlab = string('CofM=[',nnai_v_ebv_gas_perbin_centroid[0],',',$
                        nnai_v_ebv_gas_perbin_centroid[1],']',$
                        format='(A0,D0.2,A0,D0.2,A0)')
        al_legend,[pklab,comlab],col=['BLU3','RED3'],psym=[9,16],symsize=[2,2],$
                  /left,/top,box=0,spacing=1.5,thick=[4d,4d]
        cgtext,'Unweighted',-0.2,11.2,align=0.5
;        result = gauss2dfit(scaledDensity_perbin,a,xcon,ycon,/tilt)
;        cgoplot,[a[4]],[a[5]],psym=16,symsize=1,color='Orange'
;        cgoplot,[a[4]+a[2]*cos(a[6]),a[4]-a[2]*cos(a[6])],$
;                [a[5]-a[2]*sin(a[6]),a[5]+a[2]*sin(a[6])],$
;                color='Orange',thick=4
;        cgoplot,[a[4]-a[3]*sin(a[6]),a[4]+a[3]*sin(a[6])],$
;                [a[5]-a[3]*cos(a[6]),a[5]+a[3]*cos(a[6])],$
;                color='Orange',thick=4
     endif else begin
        cgimage,rebin(scaledDensity_perarea,(dsize[1]-1)*samplefac2,$
                      (dsize[2]-1)*samplefac2,/sample),$
                xran=xran,yran=yran,$
                Position=[0.55, 0.125, 0.95, 0.95], /Axes, /noerase, /keep

        tvlct,[[27],[158],[119]],100
        tvlct,[[217],[95],[2]],101
        tvlct,[[117],[112],[179]],102
        tvlct,[[231],[41],[138]],103
        tvlct,[[0],[0],[0]],104
        tvlct,[[102],[166],[30]],105
        ;  Hobbs 1974
        ;  N(NaI) = 1.7x10^14 E(B-V)^1.8
        ebv_hobbs = dindgen(26)/25d*0.5d
        logebv_hobbs = alog10(ebv_hobbs)
        logebv_hobbs_extrap = logebv_hobbs+1d
        nnai_hobbs = 14d + alog10(1.7) + 1.8d*logebv_hobbs
        nnai_hobbs_extrap = 14d + alog10(1.7) + 1.8d*logebv_hobbs_extrap
        cgoplot,logebv_hobbs,nnai_hobbs,color=105,thick=4
        cgoplot,logebv_hobbs,nnai_hobbs+0.3d,color=105,thick=2
        cgoplot,logebv_hobbs,nnai_hobbs-0.3d,color=105,thick=2
        cgoplot,logebv_hobbs_extrap,nnai_hobbs_extrap,color=105,thick=4,linesty=1
        cgoplot,logebv_hobbs_extrap,nnai_hobbs_extrap+0.3d,color=105,thick=2,linesty=1
        cgoplot,logebv_hobbs_extrap,nnai_hobbs_extrap-0.3d,color=105,thick=2,linesty=1

        ;  Phillips et al. 2013: fit to Sembach+93 + their own MW data
        ebv_sembach_phillips = dindgen(41)/40d*0.4d
        logebv_sembach_phillips = alog10(ebv_sembach_phillips)
        logebv_sembach_phillips_extrap = logebv_sembach_phillips + 1d
        nnai_sembach_phillips = 13.180d + 1.125d*(logebv_sembach_phillips+alog10(3d))
        nnai_sembach_phillips_extrap = 13.180d + 1.125d*(logebv_sembach_phillips_extrap+alog10(3d))
        cgoplot, logebv_sembach_phillips, nnai_sembach_phillips,color='Black',thick=4
        cgoplot, logebv_sembach_phillips, nnai_sembach_phillips+0.26d,color='Black',thick=2
        cgoplot, logebv_sembach_phillips, nnai_sembach_phillips-0.26d,color='Black',thick=2
        cgoplot, logebv_sembach_phillips_extrap, nnai_sembach_phillips_extrap,color='Black',thick=4,linesty=1
        cgoplot, logebv_sembach_phillips_extrap, nnai_sembach_phillips_extrap+0.26d,color='Black',thick=2,linesty=1
        cgoplot, logebv_sembach_phillips_extrap, nnai_sembach_phillips_extrap-0.26d,color='Black',thick=2,linesty=1
        ;   cgoplot, logebv_sembach_phillips, nnai_sembach_phillips+1d,color='Black',thick=8,linesty=2

        cgoplot,[-1.7,0.3],[13.05,15.55],color=101,thick=8,linesty=1
        cgoplot,[-1.7,0.3],[10.6,13.35],color=101,thick=8,linesty=1

        cgoplot,nnai_v_ebv_gas_perbin_x,nnai_v_ebv_gas_perbin_y,psym=16,symsize=0.2
        cgcontour,normDensity_perarea,xcon,ycon,/noerase,/over,c_color='Red',$
                  levels=[0.2,0.4,0.8],c_thick=4
        maxbin = max(normDensity_perarea,pixmax1d)
        pixmax = array_indices(normDensity_perarea,pixmax1d)
        nnai_v_ebv_gas_perarea_max = dblarr(2)
        nnai_v_ebv_gas_perarea_max[0] = $
          interpol(xcon,dindgen(dsize[1]-1),[pixmax[0]])
        nnai_v_ebv_gas_perarea_max[1] = $
          interpol(ycon,dindgen(dsize[2]-1),[pixmax[1]])
        cgoplot,[nnai_v_ebv_gas_perarea_max[0]],$
                [nnai_v_ebv_gas_perarea_max[1]],psym=9,symsize=2,color='BLU3',$
                err_xlow=binsizex/2d,err_xhi=binsizex/2d,$
                err_ylow=binsizey/2d,err_yhi=binsizey/2d,$
                err_thick=4d,thick=4d
        pixcent = centroid(normDensity_perarea)
        nnai_v_ebv_gas_perarea_centroid = dblarr(2)
        nnai_v_ebv_gas_perarea_centroid[0] = $
          interpol(xcon,dindgen(dsize[1]-1),[pixcent[0]])
        nnai_v_ebv_gas_perarea_centroid[1] = $
          interpol(ycon,dindgen(dsize[2]-1),[pixcent[1]])
        cgoplot,[nnai_v_ebv_gas_perarea_centroid[0]],$
                [nnai_v_ebv_gas_perarea_centroid[1]],$
                psym=16,symsize=2,color='RED3',$
                err_xlow=binsizex/2d,err_xhi=binsizex/2d,$
                err_ylow=binsizey/2d,err_yhi=binsizey/2d,$
                err_thick=4d
        pklab = string('Peak=[',nnai_v_ebv_gas_perarea_max[0],',',$
                       nnai_v_ebv_gas_perarea_max[1],']',$
                       format='(A0,D0.1,A0,D0.1,A0)')
        comlab = string('CofM=[',nnai_v_ebv_gas_perarea_centroid[0],',',$
                        nnai_v_ebv_gas_perarea_centroid[1],']',$
                        format='(A0,D0.2,A0,D0.2,A0)')
        al_legend,[pklab,comlab],col=['BLU3','RED3'],psym=[9,16],symsize=[2,2],$
                  /left,/top,box=0,spacing=1.5,thick=[4d,4d]
        cgtext,'Area-weighted',-0.2,11.2,align=0.5
;        result = gauss2dfit(scaledDensity_perarea,a,xcon,ycon,/tilt)
;        cgoplot,[a[4]],[a[5]],psym=16,symsize=1,color='Orange'
;        cgoplot,[a[4]+a[2]*cos(a[6]),a[4]-a[2]*cos(a[6])],$
;                [a[5]-a[2]*sin(a[6]),a[5]+a[2]*sin(a[6])],$
;                color='Orange',thick=4
;        cgoplot,[a[4]-a[3]*sin(a[6]),a[4]+a[3]*sin(a[6])],$
;                [a[5]-a[3]*cos(a[6]),a[5]+a[3]*cos(a[6])],$
;                color='Orange',thick=4
     endelse
  endfor
  
  cgps_close

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; NaD vel vs. stellar vel
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
     
   !P.thick=1
     
   nx = 4
   ny = 2

   cgps_open,plotdir+'s7nad_NaDvel_v_stelvel.eps',charsize=0.5,$
             /encap,/nomatch,/inches,xs=plotquantum*nx*0.95d,ys=plotquantum*ny,/qui
;  units are pixels per tenths of an inch
   pos = cglayout([nx,ny],oxmar=[6,0],oymar=[4,2],ixmar=[-2,-1],iymar=[0,0],$
                  unit=!D.X_PX_CM*2.54d/10d,xgap=4d,ygap=5d)
   noerase=0b
   xran = [-250,250]
   yran = xran
   cbdivinit=100d
   for i=0,7 do begin
      xdr = mapdir+gals[i]+'/vormerge/'+gals[i]+'_plots.xdr'
      if file_test(xdr) then begin
         restore,file=xdr
   
         if i gt 0 then noerase=1b
         cgplot,[0],/xsty,/ysty,xran=xran,yran=yran,aspect=1d,$
                pos=pos[*,i],noerase=noerase,title=plotstr.name,/nodata


;        Emission
         map = plotstr.nademvel[*,*,0]
         maperr = plotstr.errnademvel[*,*,0]
         igd_thiscomp = where(map ne bad,ctgd_thiscomp)
         ibd_thiscomp = where(map eq bad,ctbd_thiscomp)
         if ctgd_thiscomp gt 0 then $
            igd = cgsetintersection(plotstr.igd_nadem_fitweq,igd_thiscomp) $
         else igd=plotstr.igd_nadem_fitweq
         if ctbd_thiscomp gt 0 then $
            ibd = cgsetunion(plotstr.ibd_nadem_fitweq,ibd_thiscomp) $
         else ibd = plotstr.ibd_nadem_fitweq
         map[ibd] = bad
         maperr[ibd] = 0d
         cgoplot,plotstr.stel_vel,map,psym=16,symsize=0.5d,color='Grey',$
                 err_xlo=plotstr.stel_errvel[*,*,0],err_xhi=plotstr.stel_errvel[*,*,1],$
                 err_ylo=maperr[*,*,0],err_yhi=maperr[*,*,0],/err_clip,$
                 err_width=0d

;        Absorption
         if plotstr.plotabscvdf then begin
            map = plotstr.nadabscvdfvals['vpk']
            maperr = plotstr.errnadabsvel[*,*,0,*]*0d
         endif else begin
            map = plotstr.nadabsvel[*,*,0]
            maperr = plotstr.errnadabsvel[*,*,0,*]
         endelse
         igd_thiscomp = where(map ne bad,ctgd_thiscomp)
         ibd_thiscomp = where(map eq bad,ctbd_thiscomp)
         if ctgd_thiscomp gt 0 then $
            igd = cgsetintersection(plotstr.igd_nadabs_fitweq,igd_thiscomp) $
         else igd=plotstr.igd_nadabs_fitweq
         if ctbd_thiscomp gt 0 then $
            ibd = cgsetunion(plotstr.ibd_nadabs_fitweq,ibd_thiscomp) $
         else ibd = plotstr.ibd_nadabs_fitweq
         map[ibd] = bad
         maperr[ibd] = 0d
         cgoplot,plotstr.stel_vel,map,psym=16,symsize=0.5d,$
                err_xlo=plotstr.stel_errvel[*,*,0],err_xhi=plotstr.stel_errvel[*,*,1],$
                err_ylo=maperr[*,*,0,0],err_yhi=maperr[*,*,0,1],/err_clip,$
                err_width=0d
         cgoplot,xran,xran
         cgoplot,yran,yran
         cgoplot,xran,[0,0],/linesty
         cgoplot,[0,0],yran,/linesty
      
         if gals[i] eq 'eso339' then begin
            legtext = ['v$\downpeak$ (abs)','v$\down50$ (em)']
            legpsym = [16,16]
            legcol = ['Black','Grey']
            al_legend,legtext,/bottom,/right,psym=legpsym,colors=legcol,$
                      symsize=0.75d,background_col='White',chars=1
         endif

      endif
   endfor

   xtit='Stellar velocity v'+substar+' (km/s)'
   ytit='NaID velocity (km/s)'
   cgtext,xtit,0.5d,0.02d,align=0.5,/norm,chars=1
   cgtext,ytit,0.03d,0.5d,align=0.5,orient=90d,/norm,chars=1

   cgps_close

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; NaD vel -- abs minus em
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
     
   velsigthresh=2d
     
   nx = 5
   ny = 1

   cgps_open,plotdir+'s7nad_NaDvel_abs_minus_em_v_stelvel.eps',charsize=1,$
             /encap,/nomatch,/inches,xs=plotquantum*nx*0.95d,ys=plotquantum*ny,/qui
;  units are pixels per tenths of an inch
   pos = cglayout([nx,ny],oxmar=[7,0],oymar=[4,2],ixmar=[-2,-1],$
                  unit=!D.X_PX_CM*2.54d/10d,xgap=4d,ygap=5d)
   noerase=0b
   xran = [-299,299]
   yran = xran
   cbdivinit=100d
   gals_of = ['ngc1266','ngc1808','eso500','ngc5728','ic1481']
   for i=0,n_elements(gals_of)-1 do begin
      xdr = mapdir+gals_of[i]+'/vormerge/'+gals_of[i]+'_plots.xdr'
      if file_test(xdr) then begin
         restore,file=xdr
;        Absorption
         if plotstr.plotabscvdf then begin
            mapabs = plotstr.nadabscvdfvals['vpk']
            maperrabs = plotstr.errnadabsvel[*,*,0,*]*0d
         endif else begin
            mapabs = plotstr.nadabsvel[*,*,0]
            maperrabs = plotstr.errnadabsvel[*,*,0,*]
         endelse
         igd_thiscomp = where(mapabs ne bad,ctgd_thiscomp)
         ibd_thiscomp = where(mapabs eq bad,ctbd_thiscomp)
         if ctgd_thiscomp gt 0 then $
            igdabs = cgsetintersection(plotstr.igd_nadabs_fitweq,igd_thiscomp) $
         else igdabs=plotstr.igd_nadabs_fitweq
         if ctbd_thiscomp gt 0 then $
            ibd = cgsetunion(plotstr.ibd_nadabs_fitweq,ibd_thiscomp) $
         else ibd = plotstr.ibd_nadabs_fitweq
         mapabs[ibd] = bad
         maperrabs[ibd] = 0d

         mapem = plotstr.nademvel[*,*,0]
         maperrem = plotstr.errnademvel[*,*,0]
         igd_thiscomp = where(mapem ne bad,ctgd_thiscomp)
         ibd_thiscomp = where(mapem eq bad,ctbd_thiscomp)
         if ctgd_thiscomp gt 0 then $
            igdem = cgsetintersection(plotstr.igd_nadem_fitweq,igd_thiscomp) $
         else igdem=plotstr.igd_nadem_fitweq
         if ctbd_thiscomp gt 0 then $
            ibd = cgsetunion(plotstr.ibd_nadem_fitweq,ibd_thiscomp) $
         else ibd = plotstr.ibd_nadem_fitweq
         mapem[ibd] = bad
         maperrem[ibd] = 0d

         igdboth = cgsetintersection(igdabs,igdem)
         mapdiff = mapabs*0d + bad
         mapdifferrlo = mapabs*0d
         mapdifferrhi = mapabs*0d
         maperrabslo = maperrabs[*,*,0,0]
         maperrabshi = maperrabs[*,*,0,1]
         mapdiff[igdboth] = mapabs[igdboth] - mapem[igdboth]         
         mapdifferrlo[igdboth] = sqrt(maperrabslo[igdboth]^2d + maperrem[igdboth]^2d)
         mapdifferrhi[igdboth] = sqrt(maperrabshi[igdboth]^2d + maperrem[igdboth]^2d)

         stelerrlo = plotstr.stel_errvel[*,*,0]
         stelerrhi = plotstr.stel_errvel[*,*,1]

         iin = where(mapdiff ne bad $
            AND mapdiff gt 0 $
            AND (mapdiff - velsigthresh*mapdifferrlo ge 0),ctin)
         iout = where(mapdiff ne bad $
            AND mapdiff lt 0 $
            AND (mapdiff + velsigthresh*mapdifferrhi le 0),ctout)
   
         if i gt 0 then noerase=1b
         cgplot,[0],/xsty,/ysty,xran=xran,yran=yran,aspect=1d,$
                pos=pos[*,i],noerase=noerase,title=gals_of[i],/nodata,$
                charsize=0.5
         if ctin gt 0 then $
            cgoplot,plotstr.stel_vel[iin],mapdiff[iin],psym=16,symsize=0.5d,$
                    err_ylo=mapdifferrlo[iin],err_yhi=mapdifferrhi[iin],$
                    err_xlo=stelerrlo[iin],err_xhi=stelerrhi[iin],color='Red'
         if ctout gt 0 then $
            cgoplot,plotstr.stel_vel[iout],mapdiff[iout],psym=16,symsize=0.5d,$
                    err_ylo=mapdifferrlo[iout],err_yhi=mapdifferrhi[iout],$
                    err_xlo=stelerrlo[iout],err_xhi=stelerrhi[iout],color='Blue'

         cgoplot,xran,[0,0],/linesty
         cgoplot,[0,0],yran,/linesty

      
         if gals[i] eq 'ngc1266' then begin
            legtext = ['Outflow','Inflow']
            legpsym = [16,16]
            legcol = ['Blue','Red']
            al_legend,legtext,/bottom,/right,psym=legpsym,colors=legcol,$
                      symsize=0.75d,background_col='White',chars=1
         endif

      endif
   endfor

   ytit='NaD abs-em (km/s)'
   xtit='Stellar velocity (km/s)'
   cgtext,xtit,0.5d,0.02d,align=0.5,/norm,chars=1
   cgtext,ytit,0.03d,0.5d,align=0.5,orient=90d,/norm,chars=1

   cgps_close


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; NaD vel -- meta I
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;   index for total sample
   isamp = 0
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
     
   gals_of = ['ngc1266','ngc1808','eso500','ngc5728']
   ngal_of = n_elements(gals_of)
   nx = ngal_of
   ny = 5

   xsize = plotquantum*double(nx)*0.95d
   ysize = plotquantum*double(ny)
   cgps_open,plotdir+'s7nad_NaDvel_meta1.eps',charsize=1,$
             /encap,/nomatch,/inches,xs=xsize,ys=ysize,/qui
;  units are pixels per tenths of an inch
   pos = cglayout([nx,ny],oxmar=[7,0],oymar=[4,2],ixmar=[-2,-1],$
                  unit=!D.X_PX_CM*2.54d/10d,xgap=4d,ygap=5d)


   noerase=0b
   cbdivinit=100d
 
   xran = [-299,299]
   yran = xran
   for i=0,ngal_of-1 do begin
      xdr = mapdir+gals_of[i]+'/vormerge/'+gals_of[i]+'_plots.xdr'
      if file_test(xdr) then begin
         restore,file=xdr

         refvel = plotstr.stel_vel
         contours = plotstr.stel_vel_contours
         xtit='Stellar velocity v'+substar+'(km/s)'

;         line = '[NII]6583'
;         linelab = '\[NII\]6583'
;;         vtag = 'v%50'
;;         vlab = 'v50'
;         vtag = 'vpk'
;         vlab = 'vpk'
;         vplotlab = 'v$\down50$
;         emlmap = plotstr.emlvel[vtag,line]
;         ibdeml = where(emlmap eq bad OR ~ finite(emlmap),ctbdeml)
;         igdeml = where(emlmap ne bad AND emlmap ne 0 AND finite(emlmap),ctgdeml)
;         map[ibdeml] = bad   
;         refvel = emlmap
;         contours = plotstr.emlvel_contours
;         xtit='Ionized gas velocity (km/s)'

;        Absorption
         if plotstr.plotabscvdf then begin
            mapabs = plotstr.nadabscvdfvals['vpk']
            maperrabs = plotstr.errnadabsvel[*,*,0,*]*0d
         endif else begin
            mapabs = plotstr.nadabsvel[*,*,0]
            maperrabs = plotstr.errnadabsvel[*,*,0,*]
         endelse
         igd_thiscomp = where(mapabs ne bad,ctgd_thiscomp)
         ibd_thiscomp = where(mapabs eq bad,ctbd_thiscomp)
         if ctgd_thiscomp gt 0 then $
            igdabs = cgsetintersection(plotstr.igd_nadabs_fitweq,igd_thiscomp) $
         else igdabs=plotstr.igd_nadabs_fitweq
         if ctbd_thiscomp gt 0 then $
            ibd = cgsetunion(plotstr.ibd_nadabs_fitweq,ibd_thiscomp) $
         else ibd = plotstr.ibd_nadabs_fitweq
         mapabs[ibd] = bad
         maperrabs[ibd] = 0d
         mapv98abs = plotstr.nadabscvdfvals['v%98']
         mapv98abs[ibd] = bad
         mapv02abs = plotstr.nadabscvdfvals['v%02']
         mapv02abs[ibd] = bad

         mapem = plotstr.nademvel[*,*,0]
         maperrem = plotstr.errnademvel[*,*,0]
         igd_thiscomp = where(mapem ne bad,ctgd_thiscomp)
         ibd_thiscomp = where(mapem eq bad,ctbd_thiscomp)
         if ctgd_thiscomp gt 0 then $
            igdem = cgsetintersection(plotstr.igd_nadem_fitweq,igd_thiscomp) $
         else igdem=plotstr.igd_nadem_fitweq
         if ctbd_thiscomp gt 0 then $
            ibd = cgsetunion(plotstr.ibd_nadem_fitweq,ibd_thiscomp) $
         else ibd = plotstr.ibd_nadem_fitweq
         mapem[ibd] = bad
         maperrem[ibd] = 0d

         igdboth = cgsetintersection(igdabs,igdem)
         mapdiff = mapabs*0d + bad
         mapdifferrlo = mapabs*0d
         mapdifferrhi = mapabs*0d
         maperrabslo = maperrabs[*,*,0,0]
         maperrabshi = maperrabs[*,*,0,1]
         mapdiff[igdboth] = mapabs[igdboth] - mapem[igdboth]         
         mapdifferrlo[igdboth] = sqrt(maperrabslo[igdboth]^2d + maperrem[igdboth]^2d)
         mapdifferrhi[igdboth] = sqrt(maperrabshi[igdboth]^2d + maperrem[igdboth]^2d)

         igddiff = [cgsetdifference(igdabs,igdem),cgsetdifference(igdem,igdabs)]         

;        outflow / inflow from differences between em/abs only
         iin = where(mapdiff ne bad $
            AND mapdiff gt 0 $
            AND (mapdiff - velsigthresh*mapdifferrlo ge 0),ctin)
         iout = where(mapdiff ne bad $
            AND mapdiff lt 0 $
            AND (mapdiff + velsigthresh*mapdifferrhi le 0),ctout)
;        outflow from differences with stellar velocity
         ioutabs = where((mapabs-refvel le -100 $
            AND mapabs ne bad $
            AND (mapabs + velsigthresh*maperrabshi le refvel)) $
            OR (mapv98abs-refvel le -200 $
            AND mapv98abs ne bad),ctoutabs)
         ioutem = where(mapem-refvel ge 100 $
            AND mapem ne bad $
            AND (mapem - velsigthresh*maperrem ge refvel),ctoutem)
;        inflow from differences with stellar velocity
         iinabs = where(mapabs-refvel ge 100 $
            AND mapabs ne bad $
            AND (mapabs - velsigthresh*maperrabslo ge refvel),ctinabs)
         iinem = where(mapem-refvel le -100 $
            AND mapem ne bad $
            AND (mapem + velsigthresh*maperrem le refvel),ctinem)

;        outflow/inflow from all possible criteria
         if ctout gt 0 then begin
            ioutall = iout
            if ctoutabs gt 0 then begin
               ioutall = cgsetunion(ioutall,ioutabs)
               if ctoutem gt 0 then ioutall = cgsetunion(ioutall,ioutem)
            endif else if ctoutem gt 0 then ioutall = cgsetunion(ioutall,ioutem)
            ctoutall = n_elements(ioutall)
         endif else if ctoutabs gt 0 then begin
            ioutall = ioutabs
            if ctoutem gt 0 then ioutall = cgsetunion(ioutall,ioutem)
            ctoutall = n_elements(ioutall)
         endif else if ctoutem gt 0 then begin
            ioutall = ioutem
            ctoutall = n_elements(ioutall)
         endif else ctoutall = 0
         if ctin gt 0 then begin
            iinall = iin
            if ctinabs gt 0 then begin
               iinall = cgsetunion(iinall,iinabs)
               if ctinem gt 0 then iinall = cgsetunion(iinall,iinem)
            endif else if ctinem gt 0 then iinall = cgsetunion(iinall,iinem)
            ctinall = n_elements(iinall)
         endif else if ctinabs gt 0 then begin
            iinall = iinabs
            if ctinem gt 0 then iinall = cgsetunion(iinall,iinem)
            ctinall = n_elements(iinall)
         endif else if ctinem gt 0 then begin
            iinall = iinem
            ctinall = n_elements(iinall)
         endif else ctinall = 0

;        inflow / outflow, absorption only
         ctinall_abs = 0
         if ctinall gt 0 then $
            iinall_abs = cgsetintersection(iinall,igdabs,count=ctinall_abs)
         ctoutall_abs = 0
         if ctoutall gt 0 then $
            ioutall_abs = cgsetintersection(ioutall,igdabs,count=ctoutall_abs)

         mapflow = plotstr.stel_vel
         idisk = where(plotstr.stel_vel ne bad,ctdisk)
         inodata = where(plotstr.stel_vel eq bad)
         mapflow[inodata] = 0b
         mapflow[idisk] = 1b

         if ctinabs gt 0 then mapflow[iinabs] = 2b
         if ctinem gt 0 then mapflow[iinem] = 4b
         if ctoutem gt 0 then mapflow[ioutem] = 5b

         if ctoutabs gt 0 then mapflow[ioutabs] = 3b
         if ctin gt 0 then mapflow[iin] = 6b
         if ctout gt 0 then mapflow[iout] = 7b
         sizemap = size(mapflow)
         dx = sizemap[1]
         dy = sizemap[2]

;        Compute fraction of inflow/outflow area
         diskarea_nspax[isamp] = ctdisk
         diskarea_sqkpc[isamp] = double(ctdisk)*(plotstr.kpc_per_as)^2d
         maxrad_disk[isamp] = max(plotstr.map_rkpc_ifs[idisk])
         if ctinall gt 0 then begin
            infracarea[isamp] = double(ctinall)/double(ctdisk)
            maxrad_in[isamp] = max(plotstr.map_rkpc_ifs[iinall])
         endif
         if ctoutall gt 0 then begin
            outfracarea[isamp] = double(ctoutall)/double(ctdisk)
            maxrad_out[isamp] = max(plotstr.map_rkpc_ifs[ioutall])
         endif

;        Compute averages over nnai, weqabs
         nnaigd = where(plotstr.nnai ne bad)
         weqgd = where(plotstr.weq_abs_plus_em ne bad)
         if ctinall_abs gt 0 then begin
            nnaiingd = cgsetintersection(nnaigd,iinall_abs,count=nnaiinct)
            if nnaiinct gt 0 then begin
               nnaiall_in = [nnaiall_in,plotstr.nnai[nnaiingd]]
               nnaiavg_in[isamp] = mean(plotstr.nnai[nnaiingd])
               nnaisdev_in[isamp] = stddev(plotstr.nnai[nnaiingd])
            endif
         endif
         if ctinall gt 0 then begin
            weqall_in = [weqall_in,plotstr.weq_abs_plus_em[iinall]]
            weqavg_in[isamp] = mean(plotstr.weq_abs_plus_em[iinall])
            weqsdev_in[isamp] = stddev(plotstr.weq_abs_plus_em[iinall])
         endif
         if ctoutall_abs gt 0 then begin
            nnaioutgd = cgsetintersection(nnaigd,ioutall_abs,count=nnaioutct)
            if nnaioutct gt 0 then begin
               nnaiall_out = [nnaiall_out,plotstr.nnai[nnaioutgd]]
               nnaiavg_out[isamp] = mean(plotstr.nnai[nnaioutgd])
               nnaisdev_out[isamp] = stddev(plotstr.nnai[nnaioutgd])
            endif
         endif
         if ctoutall gt 0 then begin
            weqall_out = [weqall_out,plotstr.weq_abs_plus_em[ioutall]]
            weqavg_out[isamp] = mean(plotstr.weq_abs_plus_em[ioutall])
            weqsdev_out[isamp] = stddev(plotstr.weq_abs_plus_em[ioutall])
         endif
         
         mapscl = cgimgscl(rebin(mapflow,dx*samplefac,dy*samplefac,/sample),$
                                 minval=0,maxval=max(mapflow),ncolors=max(mapflow)+1)
   
         if i gt 0 then noerase=1b

;         colinabs = 'RED3'
;         coloutabs = 'BLU5'
;         colinem = 'RED3'
;         coloutem = 'BLU3'
;         colin = 'RED5'
;         colout = 'BLU5'
         colinabs = 'RED3'
         coloutabs = 'BLU3'
         colinem = 'RED3'
         coloutem = 'BLU3'
         colin = 'RED3'
         colout = 'BLU3'

;        Abs - em
         cgplot,[0],/xsty,/ysty,xran=xran,yran=yran,aspect=1d,$
                pos=pos[*,i],noerase=noerase,title=plotstr.name,/nodata,$
                charsize=0.5
         if ctin gt 0 then $
            cgoplot,refvel[iin],mapdiff[iin],psym=16,symsize=0.5d,$
                    err_ylo=mapdifferrlo[iin],err_yhi=mapdifferrhi[iin],$
                    color=colin,err_width=0d
         if ctout gt 0 then $
            cgoplot,refvel[iout],mapdiff[iout],psym=16,symsize=0.5d,$
                    err_ylo=mapdifferrlo[iout],err_yhi=mapdifferrhi[iout],$
                    color=colout,err_width=0d

         cgoplot,xran,[0,0],/linesty
         cgoplot,[0,0],yran,/linesty
         if gals[i] eq 'ngc1266' then begin
            legtext = ['Out (P Cygni)','In (Inverse P Cygni)']
            legpsym = [16,16]
            legcol = [colout,colin]
            al_legend,legtext,/bottom,/right,psym=legpsym,colors=legcol,$
               symsize=0.75d,background_col='White',chars=1
         endif

;        Abs or em
         cgplot,[0],/xsty,/ysty,xran=xran,yran=yran,aspect=1d,$
                pos=pos[*,ngal_of+i],noerase=1b,/nodata,chars=0.5
                
         cgoplot,refvel[igddiff],mapem[igddiff],psym=1,symsize=0.5d,color='Black'
         cgoplot,refvel[iinem],mapem[iinem],psym=1,symsize=0.5d,color=colinem
         cgoplot,refvel[ioutem],mapem[ioutem],psym=1,symsize=0.5d,color=coloutem
         cgoplot,refvel[iin],mapem[iin],psym=1,symsize=0.5d,color=colin
         cgoplot,refvel[iout],mapem[iout],psym=1,symsize=0.5d,color=colout

         cgoplot,refvel[igddiff],mapabs[igddiff],psym=16,symsize=0.5d
         cgoplot,refvel[iinabs],mapabs[iinabs],psym=16,symsize=0.5d,color=colinabs
         cgoplot,refvel[ioutabs],mapabs[ioutabs],psym=16,symsize=0.5d,color=coloutabs
         cgoplot,refvel[iin],mapabs[iin],psym=16,symsize=0.5d,color=colin
         cgoplot,refvel[iout],mapabs[iout],psym=16,symsize=0.5d,color=colout

         cgoplot,xran,xran
         cgoplot,yran,yran
         cgoplot,xran,[0,0],/linesty
         cgoplot,[0,0],yran,/linesty
         if gals[i] eq 'ngc1266' then begin
            legtext = ['In (em)','In (abs)']
            legpsym = [1,16]
            legcol = [colin,colin]
            al_legend,legtext,/bottom,/right,psym=legpsym,colors=legcol,$
                      symsize=0.75d,background_col='White',chars=1
            legtext = ['Out (em)','Out (abs)']
            legpsym = [1,16]
            legcol = [colout,colout]
            al_legend,legtext,/top,/left,psym=legpsym,colors=legcol,$
                      symsize=0.75d,background_col='White',chars=1
         endif

;        OF
         yran = [-599,299]
         cgplot,[0],/xsty,/ysty,xran=xran,yran=yran,aspect=1d,$
                pos=pos[*,ngal_of*2+i],noerase=1b,/nodata,chars=0.5
         cgoplot,refvel[igddiff],mapv98abs[igddiff],psym=16,symsize=0.5d
         cgoplot,refvel[iinabs],mapv98abs[iinabs],psym=16,symsize=0.5d,color=colinabs
         cgoplot,refvel[ioutabs],mapv98abs[ioutabs],psym=16,symsize=0.5d,color=coloutabs
         cgoplot,refvel[iin],mapv98abs[iin],psym=16,symsize=0.5d,color=colin
         cgoplot,refvel[iout],mapv98abs[iout],psym=16,symsize=0.5d,color=colout
         cgoplot,xran,xran
         cgoplot,xran,[0,0],/linesty
         cgoplot,[0,0],yran,/linesty
         yran = xran

;        image
;        x size of bottom row + 1 window
         dxpos = pos[2,ngal_of*3+i]-pos[0,ngal_of*3+i]
;        y size of bottom row + 1 window
         dypos = pos[3,ngal_of*3+i]-pos[1,ngal_of*3+i]
;        x size and positions of new bottom row window minus margins
         dxposnew = dxpos * 0.9d
         pos[0,ngal_of*3+i] += dxpos * 0.05d
         pos[2,ngal_of*3+i] -= dxpos * 0.05d
;        y size and positions of new bottom row window
;        have to scale from dxposnew because x/y aspect ratio not necessarily unity
;        second factor converts to IFS aspect ratio, first factor from 
;        x to y normalized units
         dyposnew = dxposnew * xsize/ysize * double(dy)/double(dx)
         pos[3,ngal_of*3+i] -= 0.02
         pos[1,ngal_of*3+i] = pos[3,ngal_of*3+i] - dyposnew
         palette = transpose([cgcolor('White',/dec,/trip),$
                              cgcolor('BLK3',/dec,/trip),$
                              cgcolor(colinabs,/dec,/trip),$
                              cgcolor(coloutabs,/dec,/trip),$
                              cgcolor(colinem,/dec,/trip),$
                              cgcolor(coloutem,/dec,/trip),$
                              cgcolor(colin,/dec,/trip),$
                              cgcolor(colout,/dec,/trip)])
         cgimage,mapscl,pos=pos[*,ngal_of*3+i],/noerase,palette=palette,opos=truepos
         cgplot,[0],xsty=5,ysty=5,position=truepos,$
                /nodata,/noerase,xran=[0,dx],yran=[0,dy],chars=0.5
         weqmap = plotstr.weq_abs
         if gals[i] eq 'ngc1266' OR $
            gals[i] eq 'ngc1808' then weqcontours = [2,4,6,8,10]
         if gals[i] eq 'eso500' then weqcontours = [2,3,4]
         if gals[i] eq 'ngc5728' then weqcontours = [1,2,3]
         if gals[i] eq 'eso339' then weqcontours = [2,4,6]
         if gals[i] eq 'ic5063' then weqcontours = [1,2,3]
         if gals[i] eq 'ic5169' then weqcontours = [1,2,3,4]
         if gals[i] eq 'ic1481' then weqcontours = [1,2,3]
         cgcontour,weqmap,dindgen(dx)+0.5,dindgen(dy)+0.5,$
                   /overplot,c_color='White',c_linesty=0,c_thick=4,$
                   levels=weqcontours,max=100,label=0
         if contours[0] ne 0 then begin
            poscon = contours(where(contours gt 0))
            zerocon = contours(where(contours eq 0))
            negcon = contours(where(contours lt 0))
            cgcontour,refvel,dindgen(dx)+0.5,dindgen(dy)+0.5,$
                      /overplot,c_linesty=2,c_thick=2,$
                      levels=negcon,max=1000,label=0
            cgcontour,refvel,dindgen(dx)+0.5,dindgen(dy)+0.5,$
                      /overplot,c_linesty=0,c_thick=8,$
                      levels=zerocon,max=1000,label=0
            cgcontour,refvel,dindgen(dx)+0.5,dindgen(dy)+0.5,$
                      /overplot,c_linesty=0,c_thick=2,$
                      levels=poscon,max=1000,label=0

         endif
         ifsf_plotaxesnuc,plotstr.xran_kpc,plotstr.yran_kpc,$
                          plotstr.center_nuclei_kpc_x,plotstr.center_nuclei_kpc_y
;      cbpos=[pos[2,0],pos[1,0],pos[2,0]+0.02,pos[3,0]]
;      ticknames = string(dindgen(ncbdiv+1)*dzran/double(ncbdiv) - $
;         (dzran - zran[1]),format=cbform)
;      cgDCBar,ncolors=3,position=cbpos,labels=ticknames,/ver,/right,$
;         charsize=1.5,spacing=0.5
;      cgtext,0.96,0.57,textoidl('# of Absorption Components'),orient=270,/normal,$
;             align=0.5

      endif
      
      isamp++
      
   endfor

   ytit='NaID v$\downpk,abs$-v$\down50%,em$ (km/s)'
   cgtext,xtit,0.53d,0.4d,align=0.5,/norm,chars=0.85
   cgtext,ytit,0.03d,0.9d,align=0.5,orient=90d,/norm,chars=0.85
   ytit='NaID v$\downpk,abs$ or v$\down50%,em$ (km/s)'
   cgtext,ytit,0.03d,0.7d,align=0.5,orient=90d,/norm,chars=0.85
   ytit='NaID v$\down98%,abs$ (km/s)'
   cgtext,ytit,0.03d,0.5d,align=0.5,orient=90d,/norm,chars=0.85
   ytit='dy (kpc)'
   xtit='dx (kpc)'
   cgtext,xtit,0.53d,0.09d,align=0.5,/norm,chars=0.85
   cgtext,ytit,0.03d,0.22d,align=0.5,orient=90d,/norm,chars=0.85

   cgps_close

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; NaD vel -- meta II
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;  Sample index isamp defined at beginning of last plot call
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
     
   gals_of = ['eso339','ic5063','ic5169','ic1481']
   ngal_of = n_elements(gals_of)
   nx = ngal_of
   ny = 5

   xsize = plotquantum*double(nx)*0.95d
   ysize = plotquantum*double(ny)
   cgps_open,plotdir+'s7nad_NaDvel_meta2.eps',charsize=1,$
             /encap,/nomatch,/inches,xs=xsize,ys=ysize,/qui
;  units are pixels per tenths of an inch
   pos = cglayout([nx,ny],oxmar=[7,0],oymar=[4,2],ixmar=[-2,-1],$
                  unit=!D.X_PX_CM*2.54d/10d,xgap=4d,ygap=5d)


   noerase=0b
   cbdivinit=100d

   xran = [-299,299]
   yran = xran
   for i=0,ngal_of-1 do begin
      xdr = mapdir+gals_of[i]+'/vormerge/'+gals_of[i]+'_plots.xdr'
      if file_test(xdr) then begin
         restore,file=xdr

         refvel = plotstr.stel_vel
         contours = plotstr.stel_vel_contours
         xtit='Stellar velocity v'+substar+'(km/s)'

;         line = '[NII]6583'
;         linelab = '\[NII\]6583'
;;         vtag = 'v%50'
;;         vlab = 'v50'
;         vtag = 'vpk'
;         vlab = 'vpk'
;         vplotlab = 'v$\down50$
;         emlmap = plotstr.emlvel[vtag,line]
;         ibdeml = where(emlmap eq bad OR ~ finite(emlmap),ctbdeml)
;         igdeml = where(emlmap ne bad AND emlmap ne 0 AND finite(emlmap),ctgdeml)
;         map[ibdeml] = bad   
;         refvel = emlmap
;         contours = plotstr.emlvel_contours
;         xtit='Ionized gas velocity (km/s)'

;        Absorption
         if plotstr.plotabscvdf then begin
            mapabs = plotstr.nadabscvdfvals['vpk']
            maperrabs = plotstr.errnadabsvel[*,*,0,*]*0d
         endif else begin
            mapabs = plotstr.nadabsvel[*,*,0]
            maperrabs = plotstr.errnadabsvel[*,*,0,*]
         endelse
         igd_thiscomp = where(mapabs ne bad,ctgd_thiscomp)
         ibd_thiscomp = where(mapabs eq bad,ctbd_thiscomp)
         if ctgd_thiscomp gt 0 then $
            igdabs = cgsetintersection(plotstr.igd_nadabs_fitweq,igd_thiscomp) $
         else igdabs=plotstr.igd_nadabs_fitweq
         if ctbd_thiscomp gt 0 then $
            ibd = cgsetunion(plotstr.ibd_nadabs_fitweq,ibd_thiscomp) $
         else ibd = plotstr.ibd_nadabs_fitweq
         mapabs[ibd] = bad
         maperrabs[ibd] = 0d
         mapv98abs = plotstr.nadabscvdfvals['v%98']
         mapv98abs[ibd] = bad
         mapv02abs = plotstr.nadabscvdfvals['v%02']
         mapv02abs[ibd] = bad

         mapem = plotstr.nademvel[*,*,0]
         maperrem = plotstr.errnademvel[*,*,0]
         igd_thiscomp = where(mapem ne bad,ctgd_thiscomp)
         ibd_thiscomp = where(mapem eq bad,ctbd_thiscomp)
         if ctgd_thiscomp gt 0 then $
            igdem = cgsetintersection(plotstr.igd_nadem_fitweq,igd_thiscomp) $
         else igdem=plotstr.igd_nadem_fitweq
         if ctbd_thiscomp gt 0 then $
            ibd = cgsetunion(plotstr.ibd_nadem_fitweq,ibd_thiscomp) $
         else ibd = plotstr.ibd_nadem_fitweq
         mapem[ibd] = bad
         maperrem[ibd] = 0d

         igdboth = cgsetintersection(igdabs,igdem)
         mapdiff = mapabs*0d + bad
         mapdifferrlo = mapabs*0d
         mapdifferrhi = mapabs*0d
         maperrabslo = maperrabs[*,*,0,0]
         maperrabshi = maperrabs[*,*,0,1]
         mapdiff[igdboth] = mapabs[igdboth] - mapem[igdboth]         
         mapdifferrlo[igdboth] = sqrt(maperrabslo[igdboth]^2d + maperrem[igdboth]^2d)
         mapdifferrhi[igdboth] = sqrt(maperrabshi[igdboth]^2d + maperrem[igdboth]^2d)

         igddiff = [cgsetdifference(igdabs,igdem),cgsetdifference(igdem,igdabs)]         

;        outflow / inflow from differences between em/abs only
         iin = where(mapdiff ne bad $
            AND mapdiff gt 0 $
            AND (mapdiff - velsigthresh*mapdifferrlo ge 0),ctin)
         iout = where(mapdiff ne bad $
            AND mapdiff lt 0 $
            AND (mapdiff + velsigthresh*mapdifferrhi le 0),ctout)
;        outflow from differences with stellar velocity
         ioutabs = where((mapabs-refvel le -100 $
            AND mapabs ne bad $
            AND (mapabs + velsigthresh*maperrabshi le refvel)) $
            OR (mapv98abs-refvel le -200 $
            AND mapv98abs ne bad),ctoutabs)
         ioutem = where(mapem-refvel ge 100 $
            AND mapem ne bad $
            AND (mapem - velsigthresh*maperrem ge refvel),ctoutem)
;        inflow from differences with stellar velocity
         iinabs = where(mapabs-refvel ge 100 $
            AND mapabs ne bad $
            AND (mapabs - velsigthresh*maperrabslo ge refvel),ctinabs)
         iinem = where(mapem-refvel le -100 $
            AND mapem ne bad $
            AND (mapem + velsigthresh*maperrem le refvel),ctinem)

;        outflow/inflow from all possible criteria
         if ctout gt 0 then begin
            ioutall = iout
            if ctoutabs gt 0 then begin
               ioutall = cgsetunion(ioutall,ioutabs)
               if ctoutem gt 0 then ioutall = cgsetunion(ioutall,ioutem)
            endif else if ctoutem gt 0 then ioutall = cgsetunion(ioutall,ioutem)
            ctoutall = n_elements(ioutall)
         endif else if ctoutabs gt 0 then begin
            ioutall = ioutabs
            if ctoutem gt 0 then ioutall = cgsetunion(ioutall,ioutem)
            ctoutall = n_elements(ioutall)
         endif else if ctoutem gt 0 then begin
            ioutall = ioutem
            ctoutall = n_elements(ioutall)
         endif else ctoutall = 0
         if ctin gt 0 then begin
            iinall = iin
            if ctinabs gt 0 then begin
               iinall = cgsetunion(iinall,iinabs)
               if ctinem gt 0 then iinall = cgsetunion(iinall,iinem)
            endif else if ctinem gt 0 then iinall = cgsetunion(iinall,iinem)
            ctinall = n_elements(iinall)
         endif else if ctinabs gt 0 then begin
            iinall = iinabs
            if ctinem gt 0 then iinall = cgsetunion(iinall,iinem)
            ctinall = n_elements(iinall)
         endif else if ctinem gt 0 then begin
            iinall = iinem
            ctinall = n_elements(iinall)
         endif else ctinall = 0

;        inflow / outflow, absorption only
         ctinall_abs = 0
         if ctinall gt 0 then $
            iinall_abs = cgsetintersection(iinall,igdabs,count=ctinall_abs)
         ctoutall_abs = 0
         if ctoutall gt 0 then $
            ioutall_abs = cgsetintersection(ioutall,igdabs,count=ctoutall_abs)

         mapflow = plotstr.stel_vel
         idisk = where(plotstr.stel_vel ne bad,ctdisk)
         inodata = where(plotstr.stel_vel eq bad)
         mapflow[inodata] = 0b
         mapflow[idisk] = 1b

         if ctinabs gt 0 then mapflow[iinabs] = 2b
         if ctinem gt 0 then mapflow[iinem] = 4b
         if ctoutem gt 0 then mapflow[ioutem] = 5b

         if ctoutabs gt 0 then mapflow[ioutabs] = 3b
         if ctin gt 0 then mapflow[iin] = 6b
         if ctout gt 0 then mapflow[iout] = 7b
         sizemap = size(mapflow)
         dx = sizemap[1]
         dy = sizemap[2]
         samplefac = 10

;        Compute fraction of inflow/outflow area
         diskarea_nspax[isamp] = ctdisk
         diskarea_sqkpc[isamp] = double(ctdisk)*(plotstr.kpc_per_as)^2d
         maxrad_disk[isamp] = max(plotstr.map_rkpc_ifs[idisk])
         if ctinall gt 0 then begin
            infracarea[isamp] = double(ctinall)/double(ctdisk)
            maxrad_in[isamp] = max(plotstr.map_rkpc_ifs[iinall])
         endif
         if ctoutall gt 0 then begin
            outfracarea[isamp] = double(ctoutall)/double(ctdisk)
            maxrad_out[isamp] = max(plotstr.map_rkpc_ifs[ioutall])
         endif

;        Compute averages over nnai, weqabs
         nnaigd = where(plotstr.nnai ne bad)
         weqgd = where(plotstr.weq_abs_plus_em ne bad)
         if ctinall_abs gt 0 then begin
            nnaiingd = cgsetintersection(nnaigd,iinall_abs,count=nnaiinct)
            if nnaiinct gt 0 then begin
               nnaiall_in = [nnaiall_in,plotstr.nnai[nnaiingd]]
               nnaiavg_in[isamp] = mean(plotstr.nnai[nnaiingd])
               nnaisdev_in[isamp] = stddev(plotstr.nnai[nnaiingd])
            endif
         endif
         if ctinall gt 0 then begin
            weqall_in = [weqall_in,plotstr.weq_abs_plus_em[iinall]]
            weqavg_in[isamp] = mean(plotstr.weq_abs_plus_em[iinall])
            weqsdev_in[isamp] = stddev(plotstr.weq_abs_plus_em[iinall])
         endif
         if ctoutall_abs gt 0 then begin
            nnaioutgd = cgsetintersection(nnaigd,ioutall_abs,count=nnaioutct)
            if nnaioutct gt 0 then begin
               nnaiall_out = [nnaiall_out,plotstr.nnai[nnaioutgd]]
               nnaiavg_out[isamp] = mean(plotstr.nnai[nnaioutgd])
               nnaisdev_out[isamp] = stddev(plotstr.nnai[nnaioutgd])
            endif
         endif
         if ctoutall gt 0 then begin
            weqall_out = [weqall_out,plotstr.weq_abs_plus_em[ioutall]]
            weqavg_out[isamp] = mean(plotstr.weq_abs_plus_em[ioutall])
            weqsdev_out[isamp] = stddev(plotstr.weq_abs_plus_em[ioutall])
         endif
 
         mapscl = cgimgscl(rebin(mapflow,dx*samplefac,dy*samplefac,/sample),$
                                 minval=0,maxval=max(mapflow),ncolors=max(mapflow)+1)
   
         if i gt 0 then noerase=1b

;        Abs - em
         cgplot,[0],/xsty,/ysty,xran=xran,yran=yran,aspect=1d,$
                pos=pos[*,i],noerase=noerase,title=plotstr.name,/nodata,$
                charsize=0.5
         if ctin gt 0 then $
            cgoplot,refvel[iin],mapdiff[iin],psym=16,symsize=0.5d,$
                    err_ylo=mapdifferrlo[iin],err_yhi=mapdifferrhi[iin],$
                    color=colin
         if ctout gt 0 then $
            cgoplot,refvel[iout],mapdiff[iout],psym=16,symsize=0.5d,$
                    err_ylo=mapdifferrlo[iout],err_yhi=mapdifferrhi[iout],$
                    color=colout

         cgoplot,xran,[0,0],/linesty
         cgoplot,[0,0],yran,/linesty
         if gals[i] eq 'ngc1266' then begin
            legtext = ['Out (P Cygni)','In (Inverse P Cygni)']
            legpsym = [16,16]
            legcol = [colout,colin]
            al_legend,legtext,/bottom,/right,psym=legpsym,colors=legcol,$
               symsize=0.75d,background_col='White',chars=1
         endif

;        Abs or em
         cgplot,[0],/xsty,/ysty,xran=xran,yran=yran,aspect=1d,$
                pos=pos[*,ngal_of+i],noerase=1b,/nodata,chars=0.5
                
         cgoplot,refvel[igddiff],mapem[igddiff],psym=1,symsize=0.5d,color='Black'
         cgoplot,refvel[iinem],mapem[iinem],psym=1,symsize=0.5d,color=colinem
         cgoplot,refvel[ioutem],mapem[ioutem],psym=1,symsize=0.5d,color=coloutem
         cgoplot,refvel[iin],mapem[iin],psym=1,symsize=0.5d,color=colin
         cgoplot,refvel[iout],mapem[iout],psym=1,symsize=0.5d,color=colout

         cgoplot,refvel[igddiff],mapabs[igddiff],psym=16,symsize=0.5d
         cgoplot,refvel[iinabs],mapabs[iinabs],psym=16,symsize=0.5d,color=colinabs
         cgoplot,refvel[ioutabs],mapabs[ioutabs],psym=16,symsize=0.5d,color=coloutabs
         cgoplot,refvel[iin],mapabs[iin],psym=16,symsize=0.5d,color=colin
         cgoplot,refvel[iout],mapabs[iout],psym=16,symsize=0.5d,color=colout

         cgoplot,xran,xran
         cgoplot,yran,yran
         cgoplot,xran,[0,0],/linesty
         cgoplot,[0,0],yran,/linesty
         if gals[i] eq 'ngc1266' then begin
            legtext = ['In (em)','In (abs)']
            legpsym = [1,16]
            legcol = [colin,colin]
            al_legend,legtext,/bottom,/right,psym=legpsym,colors=legcol,$
                      symsize=0.75d,background_col='White',chars=1
            legtext = ['Out (em)','Out (abs)']
            legpsym = [1,16]
            legcol = [colout,colout]
            al_legend,legtext,/top,/left,psym=legpsym,colors=legcol,$
                      symsize=0.75d,background_col='White',chars=1
         endif

;        OF
         yran = [-599,299]
         cgplot,[0],/xsty,/ysty,xran=xran,yran=yran,aspect=1d,$
                pos=pos[*,ngal_of*2+i],noerase=1b,/nodata,chars=0.5
         cgoplot,refvel[igddiff],mapv98abs[igddiff],psym=16,symsize=0.5d
         cgoplot,refvel[iinabs],mapv98abs[iinabs],psym=16,symsize=0.5d,color=colinabs
         cgoplot,refvel[ioutabs],mapv98abs[ioutabs],psym=16,symsize=0.5d,color=coloutabs
         cgoplot,refvel[iin],mapv98abs[iin],psym=16,symsize=0.5d,color=colin
         cgoplot,refvel[iout],mapv98abs[iout],psym=16,symsize=0.5d,color=colout
         cgoplot,xran,xran
         cgoplot,xran,[0,0],/linesty
         cgoplot,[0,0],yran,/linesty
         yran = xran

;        image
;        x size of bottom row + 1 window
         dxpos = pos[2,ngal_of*3+i]-pos[0,ngal_of*3+i]
;        y size of bottom row + 1 window
         dypos = pos[3,ngal_of*3+i]-pos[1,ngal_of*3+i]
;        x size and positions of new bottom row window minus margins
         dxposnew = dxpos * 0.9d
         pos[0,ngal_of*3+i] += dxpos * 0.05d
         pos[2,ngal_of*3+i] -= dxpos * 0.05d
;        y size and positions of new bottom row window
;        have to scale from dxposnew because x/y aspect ratio not necessarily unity
;        second factor converts to IFS aspect ratio, first factor from 
;        x to y normalized units
         dyposnew = dxposnew * xsize/ysize * double(dy)/double(dx)
         pos[3,ngal_of*3+i] -= 0.02
         pos[1,ngal_of*3+i] = pos[3,ngal_of*3+i] - dyposnew
         palette = transpose([cgcolor('White',/dec,/trip),$
                              cgcolor('BLK3',/dec,/trip),$
                              cgcolor(colinabs,/dec,/trip),$
                              cgcolor(coloutabs,/dec,/trip),$
                              cgcolor(colinem,/dec,/trip),$
                              cgcolor(coloutem,/dec,/trip),$
                              cgcolor(colin,/dec,/trip),$
                              cgcolor(colout,/dec,/trip)])
         cgimage,mapscl,pos=pos[*,ngal_of*3+i],/noerase,palette=palette,opos=truepos
         cgplot,[0],xsty=5,ysty=5,position=truepos,$
                /nodata,/noerase,xran=[0,dx],yran=[0,dy],chars=0.5
         weqmap = plotstr.weq_abs
         if gals[i] eq 'ngc1266' OR $
            gals[i] eq 'ngc1808' then weqcontours = [2,4,6,8,10]
         if gals[i] eq 'eso500' then weqcontours = [2,3,4]
         if gals[i] eq 'ngc5728' then weqcontours = [1,2,3]
         if gals[i] eq 'eso339' then weqcontours = [2,4,6]
         if gals[i] eq 'ic5063' then weqcontours = [1,2,3]
         if gals[i] eq 'ic5169' then weqcontours = [1,2,3,4]
         if gals[i] eq 'ic1481' then weqcontours = [1,2,3]
         cgcontour,weqmap,dindgen(dx)+0.5,dindgen(dy)+0.5,$
                   /overplot,c_color='White',c_linesty=0,c_thick=4,$
                   levels=weqcontours,max=100 ;,label=0
         if contours[0] ne 0 then begin
            poscon = contours(where(contours gt 0))
            zerocon = contours(where(contours eq 0))
            negcon = contours(where(contours lt 0))
            cgcontour,refvel,dindgen(dx)+0.5,dindgen(dy)+0.5,$
                      /overplot,c_linesty=2,c_thick=2,$
                      levels=negcon,max=1000,label=0
            cgcontour,refvel,dindgen(dx)+0.5,dindgen(dy)+0.5,$
                      /overplot,c_linesty=0,c_thick=8,$
                      levels=zerocon,max=1000,label=0
            cgcontour,refvel,dindgen(dx)+0.5,dindgen(dy)+0.5,$
                      /overplot,c_linesty=0,c_thick=2,$
                      levels=poscon,max=1000,label=0

         endif
         ifsf_plotaxesnuc,plotstr.xran_kpc,plotstr.yran_kpc,$
                          plotstr.center_nuclei_kpc_x,plotstr.center_nuclei_kpc_y
;      cbpos=[pos[2,0],pos[1,0],pos[2,0]+0.02,pos[3,0]]
;      ticknames = string(dindgen(ncbdiv+1)*dzran/double(ncbdiv) - $
;         (dzran - zran[1]),format=cbform)
;      cgDCBar,ncolors=3,position=cbpos,labels=ticknames,/ver,/right,$
;         charsize=1.5,spacing=0.5
;      cgtext,0.96,0.57,textoidl('# of Absorption Components'),orient=270,/normal,$
;             align=0.5

      endif
      
      isamp++
      
   endfor

   ytit='NaID v$\downpk,abs$-v$\down50%,em$ (km/s)'
   cgtext,xtit,0.53d,0.4d,align=0.5,/norm,chars=0.85
   cgtext,ytit,0.03d,0.9d,align=0.5,orient=90d,/norm,chars=0.85
   ytit='NaID v$\downpk,abs$ or v$\down50%,em$ (km/s)'
   cgtext,ytit,0.03d,0.7d,align=0.5,orient=90d,/norm,chars=0.85
   ytit='NaID v$\down98%,abs$ (km/s)'
   cgtext,ytit,0.03d,0.5d,align=0.5,orient=90d,/norm,chars=0.85
   ytit='dy (kpc)'
   xtit='dx (kpc)'
   cgtext,xtit,0.53d,0.09d,align=0.5,/norm,chars=1
   cgtext,ytit,0.02d,0.25d,align=0.5,orient=90d,/norm,chars=1

   cgps_close

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; NaD vel vs. emission-line vel
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
     
   nx = 4
   ny = 2

   line = '[NII]6583'
   linelab = '\[NII\]6583'
   vtag = 'vpk'
   vlab = 'vpk'
   vplotlab = 'v$\downpeak$

   cgps_open,plotdir+'s7nad_NaDvel_v_emlvel.eps',charsize=1,$
             /encap,/nomatch,/inches,xs=plotquantum*nx*0.95d,ys=plotquantum*ny,/qui
;  units are pixels per tenths of an inch
   pos = cglayout([nx,ny],oxmar=[6,0],oymar=[4,2],ixmar=[-2,-1],iymar=[0,0],$
                  unit=!D.X_PX_CM*2.54d/10d,xgap=4d,ygap=5d)
   noerase=0b
   xran = [-250,250]
   yran = xran
   cbdivinit=100d
   for i=0,7 do begin
      xdr = mapdir+gals[i]+'/vormerge/'+gals[i]+'_plots.xdr'
      if file_test(xdr) then begin
         restore,file=xdr
   
         if i gt 0 then noerase=1b
         cgplot,[0],/xsty,/ysty,xran=xran,yran=yran,aspect=1d,$
                pos=pos[*,i],noerase=noerase,title=gals[i],/nodata,chars=0.5

         emlmap = plotstr.emlvel[vtag,line]
         ibdeml = where(emlmap eq bad OR ~ finite(emlmap),ctbdeml)
         igdeml = where(emlmap ne bad AND emlmap ne 0 AND finite(emlmap),ctgdeml)
         map[ibdeml] = bad
;        Emission
         map = plotstr.nademvel[*,*,0]
         maperr = plotstr.errnademvel[*,*,0]
         igd_thiscomp = where(map ne bad,ctgd_thiscomp)
         ibd_thiscomp = where(map eq bad,ctbd_thiscomp)
         if ctgd_thiscomp gt 0 then $
            igd = cgsetintersection(plotstr.igd_nadem_fitweq,igd_thiscomp) $
         else igd=plotstr.igd_nadem_fitweq
         if ctbd_thiscomp gt 0 then $
            ibd = cgsetunion(plotstr.ibd_nadem_fitweq,ibd_thiscomp) $
         else ibd = plotstr.ibd_nadem_fitweq
         map[ibd] = bad
         maperr[ibd] = 0d
         cgoplot,emlmap,map,psym=16,symsize=0.5d,color='Grey',$
                 err_ylo=maperr,err_yhi=maperr,/err_clip
;        Absorption
         if plotstr.plotabscvdf then begin
            map = plotstr.nadabscvdfvals['vpk']
            maperr = plotstr.errnadabsvel[*,*,0,*]*0d
         endif else begin
            map = plotstr.nadabsvel[*,*,0]
            maperr = plotstr.errnadabsvel[*,*,0,*]
         endelse
         igd_thiscomp = where(map ne bad,ctgd_thiscomp)
         ibd_thiscomp = where(map eq bad,ctbd_thiscomp)
         if ctgd_thiscomp gt 0 then $
            igd = cgsetintersection(plotstr.igd_nadabs_fitweq,igd_thiscomp) $
         else igd=plotstr.igd_nadabs_fitweq
         if ctbd_thiscomp gt 0 then $
            ibd = cgsetunion(plotstr.ibd_nadabs_fitweq,ibd_thiscomp) $
         else ibd = plotstr.ibd_nadabs_fitweq
         map[ibd] = bad
         maperr[ibd] = 0d   
         cgoplot,emlmap,map,psym=16,symsize=0.5d,$
                err_ylo=maperr[*,*,0,0],err_yhi=maperr[*,*,0,1],/err_clip

           
         cgoplot,xran,xran
         cgoplot,yran,yran
         cgoplot,xran,[0,0],/linesty
         cgoplot,[0,0],yran,/linesty

      endif
   endfor

   xtit='Ionized gas velocity (km/s)'
   ytit='NaD velocity (km/s)'
   cgtext,xtit,0.5d,0.02d,align=0.5,/norm,chars=1
   cgtext,ytit,0.03d,0.5d,align=0.5,orient=90d,/norm,chars=1

   cgps_close


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; NaD vel vs. emission-line vel -- outflows
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
     
   nx = 5
   ny = 1

   line = '[NII]6583'
   linelab = '\[NII\]6583'
   vtag = 'vpk'
   vlab = 'vpk'
   vplotlab = 'v$\downpeak$

   cgps_open,plotdir+'s7nad_NaDvel_v_emlvel_of.eps',charsize=1,$
             /encap,/nomatch,/inches,xs=plotquantum*nx*0.95d,ys=plotquantum*ny*2,/qui
;  units are pixels per tenths of an inch
   pos = cglayout([nx,ny],oxmar=[6,1],oymar=[4,2],ixmar=[0,0],$
                  unit=!D.X_PX_CM*2.54d/10d,xgap=4d,ygap=5d)
   noerase=0b
   xran = [-299,299]
   yran = [-600,300]
   cbdivinit=100d
   gals_of = ['ngc1266','ngc1808','eso500','ngc5728','ic1481']
   for i=0,n_elements(gals_of)-1 do begin
      xdr = mapdir+gals_of[i]+'/vormerge/'+gals_of[i]+'_plots.xdr'
      if file_test(xdr) then begin
         restore,file=xdr

         emlmap = plotstr.emlvel[vtag,line]
         ibdeml = where(emlmap eq bad OR ~ finite(emlmap),ctbdeml)
         igdeml = where(emlmap ne bad AND emlmap ne 0 AND finite(emlmap),ctgdeml)
         map[ibdeml] = bad
   
         if i gt 0 then noerase=1b
         cgplot,[0],/xsty,/ysty,xran=xran,yran=yran,chars=0.5,$
                pos=pos[*,i],noerase=noerase,title=gals_of[i],/nodata

;        Emission
         mapem = plotstr.nademvel[*,*,0]
         maperrem = plotstr.errnademvel[*,*,0]
         igd_thiscomp = where(mapem ne bad,ctgd_thiscomp)
         ibd_thiscomp = where(mapem eq bad,ctbd_thiscomp)
         if ctgd_thiscomp gt 0 then $
            igd = cgsetintersection(plotstr.igd_nadem_fitweq,igd_thiscomp) $
         else igd=plotstr.igd_nadem_fitweq
         if ctbd_thiscomp gt 0 then $
            ibd = cgsetunion(plotstr.ibd_nadem_fitweq,ibd_thiscomp) $
         else ibd = plotstr.ibd_nadem_fitweq
         mapem[ibd] = bad
         maperrem[ibd] = 0d

;        Absorption
         if plotstr.plotabscvdf then begin
            mapabs = plotstr.nadabscvdfvals['vpk']
            maperrabs = plotstr.errnadabsvel[*,*,0,*]*0d
         endif else begin
            mapabs = plotstr.nadabsvel[*,*,0]
            maperrabs = plotstr.errnadabsvel[*,*,0,*]
         endelse
         igd_thiscomp = where(mapabs ne bad,ctgd_thiscomp)
         ibd_thiscomp = where(mapabs eq bad,ctbd_thiscomp)
         if ctgd_thiscomp gt 0 then $
            igd = cgsetintersection(plotstr.igd_nadabs_fitweq,igd_thiscomp) $
         else igd=plotstr.igd_nadabs_fitweq
         if ctbd_thiscomp gt 0 then $
            ibd = cgsetunion(plotstr.ibd_nadabs_fitweq,ibd_thiscomp) $
         else ibd = plotstr.ibd_nadabs_fitweq
         mapabs[ibd] = bad
         maperrabs[ibd] = 0d   
         mapv98abs = plotstr.nadabscvdfvals['v%98']
         mapv98abs[ibd] = bad

         igdboth = cgsetintersection(igdabs,igdem)
         mapdiff = mapabs*0d + bad
         mapdifferrlo = mapabs*0d
         mapdifferrhi = mapabs*0d
         maperrabslo = maperrabs[*,*,0,0]
         maperrabshi = maperrabs[*,*,0,1]
         mapdiff[igdboth] = mapabs[igdboth] - mapem[igdboth]         
         mapdifferrlo[igdboth] = sqrt(maperrabslo[igdboth]^2d + maperrem[igdboth]^2d)
         mapdifferrhi[igdboth] = sqrt(maperrabshi[igdboth]^2d + maperrem[igdboth]^2d)

         iin = where(mapdiff ne bad AND mapdiff gt 0,ctin)
         iout = where(mapdiff ne bad AND mapdiff lt 0,ctout)

         cgoplot,emlmap,mapem,psym=16,symsize=0.5d,color='Grey',$
           err_ylo=maperrem,err_yhi=maperrem,/err_clip
         cgoplot,emlmap,mapv98abs,psym=16,symsize=0.5,color='Cyan'
         cgoplot,emlmap,mapabs,psym=16,symsize=0.5d,$
                err_ylo=maperrabs[*,*,0,0],err_yhi=maperrabs[*,*,0,1],/err_clip

         cgoplot,xran,xran
         cgoplot,yran,yran
         cgoplot,xran,[0,0],/linesty
         cgoplot,[0,0],yran,/linesty
      
         if gals[i] eq 'eso500' then begin
            legtext = ['v$\downpeak$ (abs)','v$\down50$ (em)','v$\down98$ (abs)']
            legpsym = [16,16,16]
            legcol = ['Black','Red','Blue']
            al_legend,legtext,/bottom,/right,psym=legpsym,colors=legcol,$
                      symsize=0.75d,background_col='White'
         endif

      endif
   endfor

   xtit='Ionized gas velocity (km/s)'
   ytit='NaD velocity (km/s)'
   cgtext,xtit,0.5d,0.02d,align=0.5,/norm,chars=1
   cgtext,ytit,0.03d,0.5d,align=0.5,orient=90d,/norm,chars=1

   cgps_close

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Plot of spectra
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

   gals_spec = ['ngc1266','ngc1266','ngc1808','ngc1808']
;   gals_name = ['NGC 1266','NGC 1266','NGC 1808','NGC 5728']
   gals_name = ['NGC 1266','NGC 1266','NGC 1808','NGC 1808']
;   spax_x = [16,17,4,1]
;   spax_y = [10,22,19,24]
   spax_x = [16,17,4,6]
   spax_y = [10,22,19,25]

   mgilines = ['MgIb5167','MgIb5173','MgIb5184']
   linelist_mgi = ifsf_linelist(mgilines,/quiet)
   mgiavg = mean((linelist_mgi.values()).toarray())

   nadlines = ['NaD1','NaD2']
   linelist_nad = ifsf_linelist(nadlines,/quiet)
   nadavg = mean((linelist_nad.values()).toarray())

   heiline = ['HeI5876']
   linelist_hei = ifsf_linelist(heiline,/quiet)

   cgps_open,plotdir+'s7nad_spectra.eps',charsize=1,$
             /encap,/nomatch,/inches,xs=7.5d,ys=10d,/qui
;  units are pixels per tenths of an inch
   pos = cglayout([2,4],oxmar=[6,1],oymar=[5,1],ixmar=[0,0],iymar=[0,0],$
                  unit=!D.X_PX_CM*2.54d/10d,xgap=4d,ygap=5d)
   noerase=0b
   yran = [0.3,1.3]
   for i=0,n_elements(gals_spec)-1 do begin
      ix = spax_x[i]-1
      iy = spax_y[i]-1
      fitxdr = specfitdir+gals_spec[i]+'/vormerge/'+gals_spec[i]+$
               '_'+string(spax_x[i],format='(I04)')+$
               '_'+string(spax_y[i],format='(I04)')+'.xdr'
      if file_test(fitxdr) then begin

         ; get fit, initialization structure, continuum cube, NaD data, and NaD fit
         restore,file=fitxdr
         initdat = call_function('ifsf_'+gals_spec[i]+'vormerge')
         restore,initdat.startempfile
         restore,initdat.outdir+initdat.label+'.cont.xdr'
         restore,initdat.outdir+initdat.label+'.nadspec.xdr'
         restore,initdat.outdir+initdat.label+'.nadfit.xdr'

         ; read data to get error
         if not tag_exist(initdat,'datext') then datext=1 else datext=initdat.datext
         if not tag_exist(initdat,'varext') then varext=2 else varext=initdat.varext
         if not tag_exist(initdat,'dqext') then dqext=3 else dqext=initdat.dqext
         header=1b
         cube = ifsf_readcube(initdat.infile,/quiet,$
            datext=datext,varext=varext,dqext=dqext)
         err = sqrt(abs(cube.var[ix,iy,*]))

         ; wavelength, continuum data, continuum fit
         wave = struct.wave / (1d + struct.zstar)
         ydat = struct.cont_dat
         ymod = struct.cont_fit
         yerr = err[struct.fitran_indx] ;*struct.ct_rchisq
         ;
         ; MgI b region
         ;
         xran = [mgiavg-30d,mgiavg+30d]
         ixranl = value_locate(wave,xran[0])
         ixranh = value_locate(wave,xran[1])
         ; Model normalization
         ynorm = max(ymod[ixranl:ixranh])
         ymodn = ymod / ynorm
         ; Data normalized as model
         ydatn = ydat / ynorm
         ; Error normalized as model
         yerrn = yerr / ynorm
         ; Convolve template with LOSVD +/- 3sigma to get range of model errors
         ; For GD05 test only
         ; templam = reform(template.lambda,n_elements(template.lambda))
         templam = template.lambda
         smsps_hi = ifsf_convsps(templam,template.flux,struct.ct_coeff,$
                                 contcube.stel_sigma[ix,iy] + contcube.stel_sigma_err[ix,iy,1])
         smsps_lo = ifsf_convsps(templam,template.flux,struct.ct_coeff,$
                                 contcube.stel_sigma[ix,iy] - contcube.stel_sigma_err[ix,iy,0])
         ; Polynomial component of model
         poly_mod = contcube.poly_mod[ix,iy,*]
         ; Resample template to data 
         smsps_hi_resamp = ifsf_interptemp(wave,templam,smsps_hi)
         smsps_lo_resamp = ifsf_interptemp(wave,templam,smsps_lo)
         redcurve = ppxf_reddening_curve(wave,struct.ct_ebv)
         smsps_hi_resamp *= redcurve
         smsps_lo_resamp *= redcurve
         ; Add in polynomial model
         smsps_hi_poly = smsps_hi_resamp + poly_mod
         smsps_lo_poly = smsps_lo_resamp + poly_mod
         ; Normalize over wavelength range of plot using max of model
         ynorm_hi = max(smsps_hi_poly[ixranl:ixranh])
         ynorm_lo = max(smsps_lo_poly[ixranl:ixranh])
         smsps_hi_polyn = smsps_hi_poly / ynorm_hi
         smsps_lo_polyn = smsps_lo_poly / ynorm_lo
         ; Use filled polygon to display
         smspsrange = [[smsps_hi_polyn],[smsps_lo_polyn]]
         smspsmax = max(smspsrange,dim=2)
         smspsmin = min(smspsrange,dim=2)

         cgplot,[0],/xsty,/ysty,xran=xran,yran=yran,chars=0.75,$
            pos=pos[*,2*i],noerase=noerase,/nodata               
         ; error
         cgcolorfill,[wave,reverse(wave)],[ydatn-yerrn,reverse(ydatn+yerrn)],$
            color='Grey',clip=[xran[0],yran[0],xran[1],yran[1]],noclip=0
         ; data
         cgoplot,wave,ydatn,color='Black',thick=1
         ; model errors
         cgcolorfill,[wave,reverse(wave)],[smspsmin,reverse(smspsmax)],$
            color='Red',clip=[xran[0],yran[0],xran[1],yran[1]],noclip=0
         ; model
         cgoplot,wave,ymodn,color='Red'
         for j=0,n_elements(mgilines)-1 do $
            cgoplot,[linelist_mgi[mgilines[j]],linelist_mgi[mgilines[j]]],$
               yran,linesty=1
         galtext = gals_name[i]+'  [x,y]=['+string(spax_x[i],format='(I0)')+$
               ','+string(spax_y[i],format='(I0)')+']'
         sigtext = '$\sigma$$\down star$ = '+$
            string(contcube.stel_sigma[ix,iy],format='(I0)')+$
            '$\up+'+string(contcube.stel_sigma_err[ix,iy,1],format='(I0)')+$
            '$$\down-'+string(contcube.stel_sigma_err[ix,iy,0],format='(I0)')+'$'
         cgtext,galtext+'  '+sigtext,$
            xran[0]+(xran[1]-xran[0])*0.05,yran[1]-(yran[1]-yran[0])*0.1,$
            /data
         ; legend
         if i eq 0 then begin
            al_legend,['Data$\+-$1$\sigma$','Stellar model$\+-$1$\sigma$'],$
               color=['Black','Red'],$
               /bottom,/right,linestyle=[0,0],linsiz=0.33,thick=[3,3]
            cgtext,'MgI b',xran[0]+(xran[1]-xran[0])*0.05,$
               yran[0]+(yran[1]-yran[0])*0.1,/data
         endif

         ;
         ; NaI D region
         ;
         noerase=1b
         xran = [nadavg-30d,nadavg+30d]
         ixranl = value_locate(wave,xran[0])
         ixranh = value_locate(wave,xran[1])
         ; normalization of model
         ynorm = max(ymod[ixranl:ixranh])
         ; divide out polynomial normalization to NaD data
         ymodn = ymod / ynorm * poly(wave,nadcube.normpars[ix,iy,*])
         ; data
         ydatn = ydat / ynorm
         yerrn = yerr / ynorm * nadcube.errcor[ix,iy]
         ; Normalization for model errors
         ; Normalize over wavelength range of plot using max of model
         ynorm_hi = max(smsps_hi_poly[ixranl:ixranh])
         ynorm_lo = max(smsps_lo_poly[ixranl:ixranh])
         smsps_hi_polyn = smsps_hi_poly / ynorm_hi * poly(wave,nadcube.normpars[ix,iy,*])
         smsps_lo_polyn = smsps_lo_poly / ynorm_lo * poly(wave,nadcube.normpars[ix,iy,*])
         ; Use filled polygon to display
         smspsrange = [[smsps_hi_polyn],[smsps_lo_polyn]]
         smspsmax = max(smspsrange,dim=2)
         smspsmin = min(smspsrange,dim=2)

         s7specres = mean(xran)/7000d/2.3538d
         param = ifsf_nadfit2param(nadfit,ix,iy)
         modflux = ifsf_nadfcn(struct.wave,param,specres=s7specres)

         cgplot,[0],/xsty,/ysty,xran=xran,yran=yran,chars=0.75,$
                pos=pos[*,2*i+1],noerase=noerase,/nodata
         ; error
         cgcolorfill,[wave,reverse(wave)],[ydatn-yerrn,reverse(ydatn+yerrn)],$
            color='Grey',clip=[xran[0],yran[0],xran[1],yran[1]],noclip=0
         ; NaD fit
         cgoplot,wave,ymodn*modflux,color='Magenta',thick=4
         ; data
         cgoplot,wave,ydatn,color='Black'
         ; stellar model errors
         cgcolorfill,[wave,reverse(wave)],[smspsmin,reverse(smspsmax)],$
            color='Red',clip=[xran[0],yran[0],xran[1],yran[1]],noclip=0
         ; stellar model
         cgoplot,wave,ymodn,color='Red'
         for j=0,n_elements(nadlines)-1 do $
            cgoplot,[linelist_nad[nadlines[j]],linelist_nad[nadlines[j]]],$
               yran,linesty=1
         cgoplot,[linelist_hei[heiline[0]],linelist_hei[heiline[0]]],$
            yran,linesty=1

         if i eq 0 then begin
            al_legend,['Line model'],color=['Magenta'],$
               /bottom,/right,linestyle=[0],linsiz=0.33,thick=3
            cgtext,'NaI D + HeI',xran[0]+(xran[1]-xran[0])*0.05,$
               yran[0]+(yran[1]-yran[0])*0.1,/data
         endif

      endif
      
   endfor
   cgtext,'Rest Wavelength ($\angstrom$)',0.45,0.01,/norm
   cgtext,'Normalized Flux',0.03,0.45,/norm,orient=90

   cgps_close


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Histogram of sigmas
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

   cgps_open,plotdir+'s7nad_hist_lsig.eps',charsize=1,$
             /encap,/nomatch,/inches,xs=5,ys=5,/qui

   xtit='$\sigma$ (km s$\up-1$)'
   cghistoplot,alog10(allsig_abs),xtit=xtit,mininput=0d,binsize=0.15d,$
      histdata=outhistsigabs,loc=outlocsigabs,/log,min_value=0.1d,$
      ytickf='(I0)',/fill,polycolor='sky blue',datacolor='pink'
   cghistoplot,alog10(allsig_em),mininput=0d,/oplot,polycolor='dodger blue',$
      nbins=n_elements(outhistsigabs),histdata=outhistsigem,$
      loc=outlocsigem,binsize=0.15d,/log,/line_fill,orient=45d,$
      datacolor='red'

   cgps_close




;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Table of inflow/outflow properties
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

   openw,lun_tmp,plotdir+'s7nad.dat',/get_lun
   printf,lun_tmp,'# Name','Rmax(disk)','Rmax(in)','Rmax(out)','A(disk)',$
          '%A(in)','%A(out)',format='(A-11,6A11)'
   printf,lun_tmp,'# -','kpc','kpc','kpc','kpc^2','-','-',format='(A-11,6A11)'
   for i=0,7 do begin
      printf,lun_tmp,gals[i],maxrad_disk[i],maxrad_in[i],maxrad_out[i],$
             diskarea_sqkpc[i],infracarea[i],outfracarea[i],$
             format='(A-11,6D11.2)'
   endfor
   printf,lun_tmp,'#'+lineofdashes,format='(A0)'
   printf,lun_tmp,'# Name','NNaIin','','NNaIout','',$
          'Weqin','','Weqout','',format='(A-11,8A8)'
   printf,lun_tmp,'# -','<>','sdev','<>','sdev',$
          '<>','sdev','<>','sdev',format='(A-11,8A8)'
   printf,lun_tmp,'# -','cm^-2','cm^-2','cm^-2','cm^-2',$
          'A','A','A','A',format='(A-11,8A8)'
   for i=0,7 do begin
      printf,lun_tmp,gals[i],$
             nnaiavg_in[i],nnaisdev_in[i],$
             nnaiavg_out[i],nnaisdev_out[i],$
             weqavg_in[i],weqsdev_in[i],$
             weqavg_out[i],weqsdev_out[i],$
             format='(A-11,8D8.2)'
   endfor
   printf,lun_tmp,'#'+lineofdashes,format='(A0)'
   printf,lun_tmp,'Abs. line sigma: mean, median, stddev (km/s)'
   printf,lun_tmp,string(mean(allsig_abs),format='(I5)'),$
      string(median(allsig_abs),format='(I5)'),$
      string(stddev(allsig_abs),format='(I5)')
   ilow1 = where(allsig_abs le 1d,ctlow1)
   ilow5 = where(allsig_abs le 5d,ctlow5)
   printf,lun_tmp,'%bins w/ sigma <= 1: ',$
      string(double(ctlow1)/double(n_elements(allsig_abs))*100d,format='(D0.3)')
   printf,lun_tmp,'%bins w/ sigma <= 5: ',$
      string(double(ctlow5)/double(n_elements(allsig_abs))*100d,format='(D0.3)')
   printf,lun_tmp,'#'+lineofdashes,format='(A0)'
   printf,lun_tmp,'Em. line sigma: mean, median, stddev (km/s)'
   printf,lun_tmp,string(mean(allsig_em),format='(I5)'),$
      string(median(allsig_em),format='(I5)'),$
      string(stddev(allsig_em),format='(I5)')
   ilow1 = where(allsig_em le 1d,ctlow1)
   ilow5 = where(allsig_em le 5d,ctlow5)
   printf,lun_tmp,'%bins w/ sigma <= 1: ',$
      string(double(ctlow1)/double(n_elements(allsig_em))*100d,format='(D0.3)')
   printf,lun_tmp,'%bins w/ sigma <= 5: ',$
      string(double(ctlow5)/double(n_elements(allsig_em))*100d,format='(D0.3)')
   printf,lun_tmp,'#'+lineofdashes,format='(A0)'
;   ; get rid of duplicates due to Voronoi tiling
;   zall = nnaiall_in+weqall_in
;   isort = bsort(zall)
;   iuniq = uniq(zall[isort])
;   isortuniq = isort[iuniq]
;   nnaiall_in_unsort = nnaiall_in[isortuniq]
;   weqall_in_unsort = weqall_in[isortuniq]
;   iresort = bsort(nnaiall_in_unsort)
;   nnaiall_in = nnaiall_in_unsort[iresort]
;   weqall_in = weqall_in_unsort[iresort]
;   ; ... end getting rid of duplicates
;   ; get rid of duplicates due to Voronoi tiling
;   zall = nnaiall_out+weqall_out
;   isort = bsort(zall)
;   iuniq = uniq(zall[isort])
;   isortuniq = isort[iuniq]
;   nnaiall_out_unsort = nnaiall_out[isortuniq]
;   weqall_out_unsort = weqall_out[isortuniq]
;   iresort = bsort(nnaiall_out_unsort)
;   nnaiall_out = nnaiall_out_unsort[iresort]
;   weqall_out = weqall_out_unsort[iresort]
;   ; ... end getting rid of duplicates
   printf,lun_tmp,'Inflow N(NaI): mean, median, stddev (log/cm^-2)'
   printf,lun_tmp,string(mean(nnaiall_in),format='(D10.2)'),$
      string(median(nnaiall_in),format='(D10.2)'),$
      string(stddev(nnaiall_in),format='(D10.2)')
   printf,lun_tmp,'Outflow N(NaI): mean, median, stddev (log/cm^-2)'
   printf,lun_tmp,string(mean(nnaiall_out),format='(D10.2)'),$
      string(median(nnaiall_out),format='(D10.2)'),$
      string(stddev(nnaiall_out),format='(D10.2)')
   printf,lun_tmp,'Inflow total Weq: mean, median, stddev (A)'
   printf,lun_tmp,string(mean(weqall_in),format='(D10.2)'),$
      string(median(weqall_in),format='(D10.2)'),$
      string(stddev(weqall_in),format='(D10.2)')
   printf,lun_tmp,'Outflow total Weq: mean, median, stddev (A)'
   printf,lun_tmp,string(mean(weqall_out),format='(D10.2)'),$
      string(median(weqall_out),format='(D10.2)'),$
      string(stddev(weqall_out),format='(D10.2)')
   printf,lun_tmp,'#'+lineofdashes,format='(A0)'
   printf,lun_tmp,'# Fit statistics
   printf,lun_tmp,'#'+lineofdashes,format='(A0)'
   printf,lun_tmp,'# quant','r_xy','Weq0','','m_linmix','m_fitexy','dWeq_linmix',$
      'dWeq_fitexy',format='(A-8,A8,A8,A8,A10,A10,A12,A12)'
   printf,lun_tmp,'# Total weq vs stellar E(B-V)
   printf,lun_tmp,'Med',median(fit1_rxyall),median(fit1_weq0all),'',$
      median(fit1_mall_linmix),median(fit1_mall_fitexy),$
      median(fit1_dweqall_linmix),median(fit1_dweqall_fitexy),$
      format='(A-8,D8.2,D8.2,A8,D10.2,D10.2,D12.2,D12.2)'
   printf,lun_tmp,'StdDEv',stddev(fit1_rxyall),stddev(fit1_weq0all),'',$
      stddev(fit1_mall_linmix),stddev(fit1_mall_fitexy),$
      stddev(fit1_dweqall_linmix),stddev(fit1_dweqall_fitexy),$
      format='(A-8,D8.2,D8.2,A8,D10.2,D10.2,D12.2,D12.2)'
   printf,lun_tmp,'# Total weq vs gas E(B-V)
   printf,lun_tmp,'Med',median(fit2_rxyall),median(fit2_weq0all),'',$
      median(fit2_mall_linmix),median(fit2_mall_fitexy),$
      median(fit2_dweqall_linmix),median(fit2_dweqall_fitexy),$
      format='(A-8,D8.2,D8.2,A8,D10.2,D10.2,D12.2,D12.2)'
   printf,lun_tmp,'StdDev',stddev(fit2_rxyall),stddev(fit2_weq0all),'',$
      stddev(fit2_mall_linmix),stddev(fit2_mall_fitexy),$
      stddev(fit2_dweqall_linmix),stddev(fit2_dweqall_fitexy),$
      format='(A-8,D8.2,D8.2,A8,D10.2,D10.2,D12.2,D12.2)'
   printf,lun_tmp,'# quant','r_xy','gi_lin','gi_fitexy',$
      'm_linmix','m_fitexy','dWeq_linmix',$
      'dWeq_fitexy',format='(A-8,A8,A8,A8,A10,A10,A12,A12)'
   printf,lun_tmp,'# Total weq vs color
   printf,lun_tmp,'Med',median(fit3_rxyall),median(fit3_giall_linmix),$
      median(fit3_giall_fitexy),$
      median(fit3_mall_linmix),median(fit3_mall_fitexy),$
      median(fit3_dweqall_linmix),median(fit3_dweqall_fitexy),$
      format='(A-8,D8.2,D8.2,D8.2,D10.2,D10.2,D12.2,D12.2)'
   printf,lun_tmp,'StdDev',stddev(fit3_rxyall),stddev(fit3_giall_linmix),$
      stddev(fit3_giall_fitexy),$
      stddev(fit3_mall_linmix),stddev(fit3_mall_fitexy),$
      stddev(fit3_dweqall_linmix),stddev(fit3_dweqall_fitexy),$
      format='(A-8,D8.2,D8.2,D8.2,D10.2,D10.2,D12.2,D12.2)'
   printf,lun_tmp,'#'+lineofdashes,format='(A0)'
   printf,lun_tmp,'# Correlations
   printf,lun_tmp,'#'+lineofdashes,format='(A0)'
   printf,lun_tmp,'Correlation coefficients for Weq(abs) vs. Cf:'
   printf,lun_tmp,'Median = ',median(weq_abs_cf_cc),format='(A0,D0.2)'
   printf,lun_tmp,'Min = ',min(weq_abs_cf_cc),format='(A0,D0.2)'
   printf,lun_tmp,'Max = ',max(weq_abs_cf_cc),format='(A0,D0.2)'
   printf,lun_tmp,'StdDev = ',stddev(weq_abs_cf_cc),format='(A0,D0.2)'
   printf,lun_tmp,'Correlation coefficients for Weq(abs) vs. tau:'
   printf,lun_tmp,'Median = ',median(weq_abs_tau_cc),format='(A0,D0.2)'
   printf,lun_tmp,'Min = ',min(weq_abs_tau_cc),format='(A0,D0.2)'
   printf,lun_tmp,'Max = ',max(weq_abs_tau_cc),format='(A0,D0.2)'
   printf,lun_tmp,'StdDev = ',stddev(weq_abs_tau_cc),format='(A0,D0.2)'
   printf,lun_tmp,'Correlation coefficients for Weq(tot) vs. sig(abs):'
   printf,lun_tmp,'Median = ',median(weq_tot_sigabs_cc),format='(A0,D0.2)'
   printf,lun_tmp,'Min = ',min(weq_tot_sigabs_cc),format='(A0,D0.2)'
   printf,lun_tmp,'Max = ',max(weq_tot_sigabs_cc),format='(A0,D0.2)'
   printf,lun_tmp,'StdDev = ',stddev(weq_tot_sigabs_cc),format='(A0,D0.2)'
   igde = where(weq_tot_sigem_cc ne bad)
   printf,lun_tmp,'Correlation coefficients for Weq(tot) vs. sig(em):'
   printf,lun_tmp,'Median = ',median(weq_tot_sigem_cc[igde]),format='(A0,D0.2)'
   printf,lun_tmp,'Min = ',min(weq_tot_sigem_cc[igde]),format='(A0,D0.2)'
   printf,lun_tmp,'Max = ',max(weq_tot_sigem_cc[igde]),format='(A0,D0.2)'
   printf,lun_tmp,'StdDev = ',stddev(weq_tot_sigem_cc[igde]),format='(A0,D0.2)'
   printf,lun_tmp,'#'+lineofdashes,format='(A0)'
   printf,lun_tmp,'# Stellar error statistics
   printf,lun_tmp,'#'+lineofdashes,format='(A0)'
   printf,lun_tmp,'Median/min/max of median error in vel = ',$
      median(med_stel_errvel),min(med_stel_errvel),max(med_stel_errvel),$
      format='(A0,D8.2,D8.2,D8.2)'
   printf,lun_tmp,'Median/min/max of median error in sig = ',$
      median(med_stel_errsig),min(med_stel_errsig),max(med_stel_errsig),$
      format='(A0,D8.2,D8.2,D8.2)'
   printf,lun_tmp,'Median/min/max of median error in E(B-V) = ',$
      median(med_stel_errebv),min(med_stel_errebv),max(med_stel_errebv),$
      format='(A0,D8.2,D8.2,D8.2)'
   free_lun,lun_tmp


   openw,lun_tmp,tabdir+'table_radiiareas.dat',/get_lun
   for i=0,7 do begin
     xdr = mapdir+gals[i]+'/vormerge/'+gals[i]+'_plots.xdr'
     if file_test(xdr) then begin
        restore,file=xdr
        printf,lun_tmp,plotstr.name,amp,$
               maxrad_disk[i],amp,$
               maxrad_in[i],amp,$
               maxrad_out[i],amp,$
               diskarea_sqkpc[i],amp,$
               infracarea[i]*100d,amp,$
               outfracarea[i]*100d,dslash,$
               format='(A10,A0,D6.2,A0,D6.2,A0,D6.2,A0,'+$
                      'D6.2,A0,I6.0,A0,I6.0,A0)'
      endif
   endfor
   free_lun,lun_tmp

   openw,lun_tmp,tabdir+'table_nnaiweqstats.dat',/get_lun
   for i=0,7 do begin
     xdr = mapdir+gals[i]+'/vormerge/'+gals[i]+'_plots.xdr'
     if file_test(xdr) then begin
       restore,file=xdr
       printf,lun_tmp,plotstr.name,amp,$
              nnaiavg_in[i],amp,$
              nnaisdev_in[i],amp,$
              nnaiavg_out[i],amp,$
              nnaisdev_out[i],amp,$
              weqavg_in[i],amp,$
              weqsdev_in[i],amp,$
              weqavg_out[i],amp,$
              weqsdev_out[i],dslash,$
              format='(A10,A0,D6.2,A0,D6.2,A0,D6.2,A0,'+$
                     'D6.2,A0,D6.2,A0,D6.2,A0,D6.2,A0,D6.2,A0)'
     endif
   endfor
   free_lun,lun_tmp


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; NaD: C_f vs. color
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
;   cgps_open,plotdir+'s7nad_cf_v_col.eps',charsize=1,$
;             /encap,/nomatch,/inches,xs=plotquantum*nx,ys=plotquantum*ny,/qui
;;  units are pixels per tenths of an inch
;   pos = cglayout([nx,ny],oxmar=[4,0],oymar=[4,3],ixmar=[0,0],$
;                  unit=!D.X_PX_CM*2.54d/10d,xgap=4d,ygap=5d)
;   noerase=0b
;   xran = [1.2,2.6]
;   yran = [-0.05,1.05]
;   for i=0,2 do begin
;      if i gt 0 then noerase=1b
;      xdr = mapdir+gals_nh[i]+'/vormerge/'+gals_nh[i]+'_plots.xdr'
;      restore,file=xdr
;      cmap = plotstr.cmap
;      mapabs = plotstr.cf
;;      mapabserrlo = plotstr.cf_el
;;      mapabserrhi = plotstr.cf_eh
;      igdone = plotstr.igdone_cf
;      istone = plotstr.istone_cf
;      igdtwo = plotstr.igdtwo_cf
;      isttwo = plotstr.isttwo_nh
;      cgplot,cmap[igdone],mapabs[igdone],/xsty,/ysty,xran=xran,yran=yran,$
;             color='Gray',noerase=noerase,pos=pos[*,i],$
;             psym=3,aspect=1d,err_width=0,err_color='Gray',/err_clip,$
;             err_ylow=mapabserrlo[igdone],err_yhigh=mapabserrhi[igdone],$
;             title=gals_nh[i]
;      cgoplot,cmap[istone],mapabs[istone],color='Gray',psym=3,$
;              err_width=0,err_color='Gray',/err_clip,$
;              err_ylow=mapabserrlo[istone],err_yhigh=mapabserrhi[istone]
;      if igdtwo[0] ne -1 then begin
;         cgoplot,cmap[igdtwo],mapabs[igdtwo],color='Gray',psym=3,$
;                 err_width=0,err_color='Gray',/err_clip,$
;                 err_ylow=mapabserrlo[igdtwo],err_yhigh=mapabserrhi[igdtwo]
;         cgoplot,cmap[isttwo],mapabs[isttwo],color='Gray',psym=3,$
;                 err_width=0,err_color='Gray',/err_clip,$
;                 err_ylow=mapabserrlo[isttwo],err_yhigh=mapabserrhi[isttwo]
;      endif
;      cgoplot,cmap[igdone],mapabs[igdone],color='Red',psym=9,symsize=0.75d
;      if igdtwo[0] ne -1 then $
;         cgoplot,cmap[igdtwo],mapabs[igdtwo],psym=16,symsize=0.75d
;      plotsym,2,thick=2
;      cgoplot,cmap[istone],mapabs[istone],color='Red',psym=8
;      if igdtwo[0] ne -1 then $
;         cgoplot,cmap[isttwo],mapabs[isttwo],psym=8
;      av=dindgen(101)/100d*10d
;      nh=alog10(av) + 21d + alog10(1.8d)
;      ebv_rv325 = av/3.25d
;      ebv_rv405 = av/4.05d
;      ebv_rv485 = av/4.85d
;      calz_unred,[plotstr.pivotbl,plotstr.pivotrd],[1d,1d],1d,funred
;      m = 1d/(2.5d*(alog10(funred[0])-alog10(funred[1])))
;      b = -m*plotstr.int_color
;      col_rv325 = (ebv_rv325 - b)/m
;      col_rv405 = (ebv_rv405 - b)/m
;      col_rv485 = (ebv_rv485 - b)/m
;      cgoplot,col_rv325,nh,linesty=1
;      cgoplot,col_rv405,nh
;      cgoplot,col_rv485,nh,linesty=2
;   endfor
;;   xtit='g-i',ytit='log[ N(H) / cm$\up-2$ ]',$
;   
;   cgps_close

end