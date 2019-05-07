; docformat = 'rst'
;
;+
;
; :Categories:
;    DRTOOLS
;
; :Returns:
;
; :Params:
;    initdat: in, required, type=structure
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
;      2018oct19, DSNR, created
;    
; :Copyright:
;    Copyright (C) 2018 David S. N. Rupke
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
function drt_runlinmix,x,y,xerr=xerr,yerr=yerr,yfitseed=seed,onesig=onesig,$
                       alpha=alpha,beta=beta,corr=corr,sigsqr=sigsqr

;   if ~ keyword_set(xerr) then xerr=0
;   if ~ keyword_set(yerr) then yerr=0
;   if ~ keyword_set(seed) then seed = 12345

   if ~ keyword_set(xerr) AND keyword_set(yerr) then $
      LINMIX_ERR,x,y,post,ysig=yerr,/silent $;,seed=seed
   else if ~ keyword_set(yerr) AND keyword_set(xerr) then $
      LINMIX_ERR,x,y,post,xsig=xerr,/silent $;,seed=seed
   else if ~ keyword_set(yerr) AND ~ keyword_set(xerr) then $
      LINMIX_ERR,x,y,post,/silent $ ;,seed=seed
   else $
      LINMIX_ERR,x,y,post,xsig=xerr,ysig=yerr,/silent ;,seed=seed
   npost = n_elements(post.alpha)
   nfit = n_elements(x)
   signum = fix(0.34*npost)
   twosignum = fix((1d - 0.9545d)/2d*npost)

   bestalpha = median(post.alpha)
   sortalpha = post(sort(post.alpha)).alpha
   ibestalpha = value_locate(sortalpha,bestalpha)
   sigloalpha = bestalpha - sortalpha[ibestalpha-signum]
   sighialpha = sortalpha[ibestalpha+signum] - bestalpha
   twosigloalpha = sortalpha[twosignum-1]
   twosighialpha = sortalpha[(npost-1)-twosignum]
   alpha=[twosigloalpha,bestalpha,twosighialpha]

   bestbeta = median(post.beta)
   sortbeta = post(sort(post.beta)).beta
   ibestbeta = value_locate(sortbeta,bestbeta)
   siglobeta = bestbeta - sortbeta[ibestbeta-signum]
   sighibeta = sortbeta[ibestbeta+signum] - bestbeta
   twosiglobeta = sortbeta[twosignum-1]
   twosighibeta = sortbeta[(npost-1)-twosignum]
   beta=[twosiglobeta,bestbeta,twosighibeta]

   bestcorr = median(post.corr)
   sortcorr = post(sort(post.corr)).corr
   ibestcorr = value_locate(sortcorr,bestcorr)
   siglocorr = bestcorr - sortcorr[ibestcorr-signum]
   sighicorr = sortcorr[ibestcorr+signum] - bestcorr
   twosiglocorr = sortcorr[twosignum-1]
   twosighicorr = sortcorr[(npost-1)-twosignum]
   corr=[twosiglocorr,bestcorr,twosighicorr]

   bestsigsqr = median(post.sigsqr)
   sortsigsqr = post(sort(post.sigsqr)).sigsqr
   ibestsigsqr = value_locate(sortsigsqr,bestsigsqr)
   siglosigsqr = bestsigsqr - sortsigsqr[ibestsigsqr-signum]
   sighisigsqr = sortsigsqr[ibestsigsqr+signum] - bestsigsqr
   twosiglosigsqr = sortsigsqr[twosignum-1]
   twosighisigsqr = sortsigsqr[(npost-1)-twosignum]
   sigsqr=[twosiglosigsqr,bestsigsqr,twosighisigsqr]

   fitvals = median(post.alpha) + median(post.beta)*x
   modpts = rebin(post.alpha,npost,nfit) + $
            rebin(post.beta,npost,nfit)*$
            rebin(transpose(x),npost,nfit)
   yloenv = dblarr(nfit)
   yhienv = dblarr(nfit)
   for i=0,nfit-1 do begin
      sortvals = modpts(sort(modpts(*,i)),i)
      if keyword_set(onesig) then begin
         yloenv[i] = sortvals[signum-1]
         yhienv[i] = sortvals[(npost-1)-signum]
      endif else begin
         yloenv[i] = sortvals[twosignum-1]
         yhienv[i] = sortvals[(npost-1)-twosignum]
      endelse
   endfor

   return,[[fitvals],[yloenv],[yhienv]]

end