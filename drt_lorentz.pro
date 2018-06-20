;
; History
;   13may16  DSNR  created
;

; Compute the Lorentzian function using the parameterization that's in
; MPFITPEAK.

;   Lorentzian#
;
;   Model     A[0]/(u^2 + 1)
;
;   A[0]         Peak Value
;   A[1]        Peak Centroid
;   A[2]            HWHM%
;   A[3]         + A[3]
;   A[4]         + A[4]*x
;
;   Notes: # u = (x - A[1])/A[2]
;          % Half-width at half maximum
;          * Optional depending on NTERMS

function drt_lorentz,x,a

  u = (x - a[1])/a[2]
  lor = a[0]/(u^2d + 1d)
  if n_elements(a) ge 4 then lor = lor + a[3]
  if n_elements(a) eq 5 then lor = lor + a[4]*x

  return,lor

end
