; docformat = 'rst'
;
;+
; 
; Run LINMIX_ERR with default settings and process posterior distributions to 
; extract best fit and 1sigma error envelope at each input x-value. Optionally
; return also 2sigma error envelope.
;
; :Categories:
;    DRTOOLS
;
; :Returns:
;     [5, npts] array. First column is best-fit values, second is bottom of 
;     one-sigma envelope, third is top of one-sigma envelope, fourth is 
;     bottom of two-sigma envelope, fifth is top of two-sigma envelope.
;
; :Params:
;     x: in, required, type=dblarr(N)
;     y: in, required, type=dblarr(N)
;     xerr: in, optional, type=dblarr(N)
;     yerr: in, optional, type=dblarr(N)
;    
; :Keywords:
;     alpha: out, optional, type=dblarr(5)
;        Bestfit, -1sigma, +1sigma, -2sigma, and +2sigma values for slope.
;     beta: out, optional, type=dblarr(5)
;        Same, but for intercept.
;     corr: out, optional, type=dblarr(5)
;        Same, but for correlation coefficient.
;     sigsqr: out, optional, type=dblarr(5)
;        Same, but for intrinsic scatter.
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
;      2018oct19, DSNR, created
;      2021jan28, DSNR, changed default to return
;      2021feb08, DSNR, added option for metropolitan-hastings sampler
;      2021feb09, DSNR, option to output p-value, input miniter
;    
; :Copyright:
;    Copyright (C) 2018--2021 David S. N. Rupke
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
function drt_procpostlinmix,dist,signum,twosignum,pval=pval

   npost = n_elements(dist)
   signum = fix(erf(1d/sqrt(2d))/2d*npost)
   twosignum = fix(erf(2d/sqrt(2d))/2d*npost)

   bestpar = median(dist)
   sortpar = dist(sort(dist))
   ibestpar = value_locate(sortpar,bestpar)
   siglopar = bestpar - sortpar[ibestpar-signum]
   sighipar = sortpar[ibestpar+signum] - bestpar
   twosiglopar = bestpar - sortpar[ibestpar-twosignum]
   twosighipar = sortpar[ibestpar+twosignum] - bestpar
   par=[bestpar,siglopar,sighipar,twosiglopar,twosighipar]

   ; Compute p-value if this is a correlation coefficient
   ; second element of array indicates whether or not this is an upper limit
   ; 1 = upper limit, 0 = not
   if keyword_set(pval) then begin
      ; where does 0-value appear
      izero = value_locate(sortpar,0d)
      npts_dbl = double(n_elements(sortpar))
      pvalul =  1d/npts_dbl ; default upper limit
      pval = pvalul
      ; case of positive cc
      if bestpar gt 0 then begin
         if izero gt -1 then pval = double(izero)/npts_dbl
      ; and negative cc
      endif else begin
         if izero lt n_elements(sortpar) then pval = 1d - double(izero)/npts_dbl
      endelse
      if pval eq pvalul then pval = [pvalul,1d] else pval = [pval,0d]
   endif  
     
   return,par

end


function drt_runlinmix,x,y,xerr=xerr,yerr=yerr,$
                       alpha=alpha,beta=beta,corr=corr,sigsqr=sigsqr,$
                       metro=metro,pval=pval,miniter=miniter,detected=detected,$
                       maxiter=maxiter,ngauss=ngauss,verbose=verbose

   if ~ keyword_set(metro) then metro=0b else metro=1b
   if ~ keyword_set(miniter) then miniter=5000L ; default from linmix_err
   if ~ keyword_set(ngauss) then ngauss=3 ; default from linmix_err
   if ~ keyword_set(maxiter) then miniter=100000L ; default from linmix_err
   if ~ keyword_set(detected) then detected=bytarr(n_elements(x))+1b
   if keyword_set(verbose) then silent=0b else silent=1b

   if ~ keyword_set(xerr) AND keyword_set(yerr) then $
      LINMIX_ERR,x,y,post,ysig=yerr,silent=silent,metro=metro,miniter=miniter,$
         delta=detected,maxiter=maxiter,ngauss=ngauss $
   else if ~ keyword_set(yerr) AND keyword_set(xerr) then $
      LINMIX_ERR,x,y,post,xsig=xerr,silent=silent,metro=metro,miniter=miniter,$
      delta=detected,maxiter=maxiter,ngauss=ngauss $
   else if ~ keyword_set(yerr) AND ~ keyword_set(xerr) then $
      LINMIX_ERR,x,y,post,silent=silent,metro=metro,miniter=miniter,$
      delta=detected,maxiter=maxiter,ngauss=ngauss $
   else $
      LINMIX_ERR,x,y,post,xsig=xerr,ysig=yerr,silent=silent,metro=metro,$
      miniter=miniter,delta=detected,maxiter=maxiter,ngauss=ngauss
   nfit = n_elements(x)
   npost = n_elements(post.alpha)
   signum = fix(erf(1d/sqrt(2d))/2d*npost)
   twosignum = fix(erf(2d/sqrt(2d))/2d*npost)

   alpha = drt_procpostlinmix(post.alpha)
   beta = drt_procpostlinmix(post.beta)
   pval=1b
   corr = drt_procpostlinmix(post.corr,pval=pval)
   sigsqr = drt_procpostlinmix(post.sigsqr)

   fitvals = median(post.alpha) + median(post.beta)*x
   modpts = rebin(post.alpha,npost,nfit) + $
            rebin(post.beta,npost,nfit)*$
            rebin(transpose(x),npost,nfit)
   yloenv_sig = dblarr(nfit)
   yhienv_sig = dblarr(nfit)
   yloenv_twosig = dblarr(nfit)
   yhienv_twosig = dblarr(nfit)
   for i=0,nfit-1 do begin
      sortvals = modpts(sort(modpts(*,i)),i)
      ifitval = n_elements(sortvals)/2 ; this is about the median index
      yloenv_sig[i] = sortvals[ifitval - signum]
      yhienv_sig[i] = sortvals[ifitval + signum]
      yloenv_twosig[i] = sortvals[ifitval - twosignum]
      yhienv_twosig[i] = sortvals[ifitval + twosignum]
   endfor

   return,[[fitvals],[yloenv_sig],[yhienv_sig],[yloenv_twosig],[yhienv_twosig]]

end