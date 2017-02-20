;+
;  NAME:
;    fit_rotations
;
;  PURPOSE:
;    Randomly sample clustering amplitudes from jackknife iterations
;    and fit the bias, to build a distribuiton of bias values.
;    Cross-between usual and Bayesian fitting, without relying on a
;    prior. Relies on several files from ang_cluster.pro and
;    ang_jackknife.pro, which can be names explicitly if desired.
;
;  USE:
;    fit_rotations,ell,cl,b,[optional inputs]

;  INPUTS:
;    ell - bin centers 
;    cl - cross-correlation values
;    mod_ell - model bin centers
;    mod_cl - model cls
;
;  OPTIONAL INPUTS:
;    n - number of random draws/fits (default 10000)
;    minscale - set lower limit for fitting range (default min(ell))
;    maxscale - set upper limit for fitting range (default max(ell))
;    allfile - text file with all results (5 columns, bin edge,
;                bin cent, bin width, cl, iteration; default xcorr_noise_binned.txt)
;    path - path to directory that has files (makes easier to use
;    covar - covariance matrix 
;    err - array of 1-sigma errors (used if covfile not set)
;           default names; defaults to current directory)
;    outfile - if set, will print the bias and error of each fit
;              iteration to this file
;    seed - the random seed for the random jackknife iteration draw
;           (default 818)
;
;  OUTPUTS:
;    bias - distribution of bias values
;
;  HISTORY:
;    2-12-17 - Written - MAD (Dartmouth)
;-
PRO fit_rotations,ell,cl,errs,mod_ell,mod_cl,b,n=n,minscale=minscale,maxscale=maxscale,$
                   covar=covar,$
                   outfile=outfile,seed=seed
  
  ;MAD Start timer
  st=timer()

  ;MAD Set defaults
  IF (n_elements(n) EQ 0) THEN n=1000
  IF (n_elements(n) EQ 0) THEN acc=1e-3
  IF (n_elements(seed) EQ 0) THEN seed=818
  IF (n_elements(path) EQ 0) THEN path='./'
  IF (n_elements(minscale) EQ 0) THEN minscale=min(ell)
  IF (n_elements(maxscale) EQ 0) THEN maxscale=max(ell)
  IF (n_elements(allfile) EQ 0) THEN allfile='xcorr_noise_binned.txt'
  
  ;MAD read in necessary data, convert as necessary
  readcol,allfile,blah,ellall,blah,clall,pix,format='D'

  ;MAD Find unique ell bins
  ells=ellall[rem_dup(ellall)]

  ;Initialize outputs
  b=fltarr(n)
  b_err=fltarr(n)

  ;Start n iterations
  FOR i=0,n-1 DO BEGIN
     counter,i,n
     c_l=dblarr(n_elements(ell))
     ;MAD pick a cl for each bin randomly from rotations
     FOR j=0,n_elements(ell)-1 DO BEGIN
        xx=where(ellall EQ ell[j])
        ;Have to shift to 0, widen distribution, shift back
        cls=cl[j]-clall[xx]
        indx=floor(randomu(seed)*max(pix))
        c_l[j]=cls[indx]
     ENDFOR
     ;MAD Fit random iteration
     IF (keyword_set(covar)) THEN $
        fit_crosspower,c_l,ell,minscale,maxscale,mod_ell,mod_cl,$
                       t_bias,t_biaserr,err=errs,covar=covar,/silent ELSE $
                          fit_crosspower,c_l,ell,minscale,maxscale,mod_ell,mod_cl,$
                                         t_bias,t_biaserr,err=errs,/silent
     b[i]=t_bias
     b_err[i]=t_biaserr
     
     ;MAD Hack, because readcol will die after many
     ;MAD iterations because it doesn't close luns for some reason
     close,/all
  ENDFOR

  ;MAD Write out, if needed
  IF (n_elements(outfile) NE 0) THEN BEGIN
     openw,lun,outfile,/get_lun
     printf,lun,';b       b_err',format='(A)'
     FOR i=0L,n-1 DO $
        printf,lun,strtrim(b[i],2) + '     ' + strtrim(b_err[i],2),format='(A)'
     free_lun,lun
  ENDIF
  
  et=timer(st=st,/fin)
  return
END
