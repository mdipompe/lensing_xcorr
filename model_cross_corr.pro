;+
;  NAME:
;    model_cross_corr
;  PURPOSE:
;    Generate a model DM cross-power spectrum (bias=1)
;
;  USE:
;    model_cross_corr,mod_ell,mod_cl,power_spec='power_spec.fits',dndz='dndz.txt',zarray=zarray,$
;                      omega_m=omega_m,omega_b=omega_b,omega_l=omega_l,h0=h0,$
;                      outfile='model_power.txt',paramfile='params.ini'
;
;  INPUT:
;
;  Optional Inputs:
;    power_spec - structure or string name of fits file 
;                 containing the matter power
;                 spectrum.  Needs tags pk (the power spectrum), k
;                 (wavenumber), l (angular wavenumber), z (redshift).
;                 Can be made from CAMB data using camb4idl and
;                 combine_camb.pro. If not set, will call necessary 
;                 procedures to make it.
;    paramfile - if power_spec not set, need param file for CAMB.
;                Defaults to 'default_params.ini'
;    dndz - string name of text file with dndz.  Will fit a spline
;           function and write it out to dndz_fit.txt.  You should
;           check this!!  Defaults to dndz.txt
;    zarray - array of z values power spectrum and dndz are calculated
;             for. If not set, defaults to 0.01 to 4.0 in steps of 0.01
;    omega_m - Omega_matter, defaults to 0.273.  Be sure it is
;              consistent with what you used when you calculated the
;              power spectrum if you supply one!
;    omega_l - Omega_lambda, defaults to 0.727.  See above...
;    h0 - little h (H0/100), defaults to 0.702.  See above...
;    omega_b - omega_baryon, defaults to 0.046
;    outfile - if supplied, writes model power out to text file
;    plotfile - if supplied, makes a plot of the model
;    bz - coefficients of b(z) model, of the form b(z) = bz[0] + bz[1](1+z)^2
;
;  KEYWORDS:
;    
;  OUTPUT:
;    mod_ell - the ell values of the model
;    mod_cl - the model cross power
;    dndz_fit.txt - a file containing the normalized dndz fit as a check
;
;  NOTES:
;    Make sure your cosmology agrees with how you calculated the power
;    spectrum, if you supply your own. Also check the dndz_fit.txt file to make sure nothing
;    went wrong!
;
;  HISTORY:
;    11-7-14 - Written - MAD (UWyo)
;    8-21-15 - If power spec not supplied, calls CAMB4IDL to get it -
;              MAD (UWyo)
;    11-6-15 - Added bz option - MAD (Dartmouth)
;-
PRO model_cross_corr,mod_ell,mod_cl,power_spec=power_spec,dndz=dndz,zarray=zarray,$
                     paramfile=paramfile,$
                     omega_m=omega_m,omega_b=omega_b,omega_l=omega_l,h0=h0,$
                     outfile=outfile,plotfile=plotfile,bz=bz

;MAD If output file already exists, don't run just read it in
IF ~keyword_set(outfile) THEN check='' ELSE check=file_search(outfile)
IF (check NE '') THEN BEGIN
   print,'MODEL_CROSS_CORR: Output file already exists, reading in and returning...'
   readcol,outfile,mod_ell,mod_cl,format='F,D'
   return
ENDIF

;MAD Set some defaults
IF ~keyword_set(dndz) THEN dndz='dndz.txt'
IF ~keyword_set(zarray) THEN zarray=(findgen(400)/100.)+0.01
IF ~keyword_set(paramfile) THEN paramfile='default_params.ini'

c=2.99792458e5
IF ~keyword_set(h0) THEN h0=0.702
IF ~keyword_set(omega_m) THEN omega_m=0.275
IF ~keyword_set(omega_b) THEN omega_b=0.046
IF ~keyword_set(omega_l) THEN omega_l=0.725

;MAD Convert z array to comoving distance array
chi=fltarr(n_elements(zarray))
FOR i=0L,n_elements(zarray)-1 DO BEGIN
   d=cosmocalc(zarray[i],h=h0,om=omega_m,lambda=omega_l)
   chi[i]=d.d_c
ENDFOR
;chi=cosmo_gen('c',omega_m,omega_l,h0,zarray) 


IF ~keyword_set(power_spec) THEN BEGIN
   root='camb'
   revz=reverse(zarray)
   numloop=floor((n_elements(revz)/150))
   IF (numloop NE 0) THEN BEGIN
      numrem=(n_elements(revz) MOD 150)
      FOR i=0L,numloop-1 DO BEGIN
         indx=indgen(150)+(i*150)
         matter_power_spec,paramfile,zarray[indx],h0=h0,omega_b=omega_b,omega_dm=omega_m-omega_b,omega_l=omega_l
      ENDFOR
      IF (numrem NE 0) THEN BEGIN
         matter_power_spec,paramfile,zarray[(n_elements(zarray)-numrem):n_elements(zarray)-1],$
                           h0=h0,omega_b=omega_b,omega_dm=omega_m-omega_b,omega_l=omega_l
      ENDIF
   ENDIF ELSE BEGIN
      matter_power_spec,paramfile,zarray,h0=h0,omega_b=omega_b,omega_dm=omega_m-omega_b,omega_l=omega_l
   ENDELSE
   combine_camb,'./',zarray,pspec,outfile='power_spec_camb.fits',maxell=3000.,chiin=chi*h0
ENDIF ELSE BEGIN
;MAD Read in power spectrum if given file name instead of structure
   check=size(power_spec)
   IF ((check[0] EQ 0) AND (check[1] EQ 7)) THEN $
      pspec=mrdfits(power_spec,1) ELSE $
         pspec=power_spec
ENDELSE

;Put h factors in power spectrum
delsq=pspec.pk*(1./(2.*!dpi^2.))*pspec.k^3.
pkk=pspec.k*h0
pkl=pspec.l
pkz=pspec.z
pk=delsq*(2.*!dpi^2.)/(pkk^3.)
undefine,pspec


;MAD Read in redshift distributions
readcol,dndz,zdist,format='D'

;MAD Get comoving distance to CMB
;d_cmb=cosmo_gen('c',omega_m,omega_l,h0,1100.)
d_cmb=cosmocalc(1100.,om=omega_m,lambda=omega_l,h=h0)
d_cmb=d_cmb.d_c

;MAD Get array of k, l values
kvals=pkk[where(pkz EQ min(pkz))]
lvals=pkl[where(pkz EQ min(pkz))]

;MAD Fit z values to get dndz
fit_dndz,zdist,zarray,dndz

;MAD Do the integral part of C_l
c_l=fltarr(n_elements(lvals))
FOR i=0L,n_elements(c_l)-1 DO Begin
   xx=where(pkl EQ lvals[i])
   IF ~keyword_set(bz) THEN fz=((d_cmb[0]-chi)/(d_cmb[0]))*((1.+zarray)/chi)*dndz*pk[xx] ELSE $
      fz=((d_cmb[0]-chi)/(d_cmb[0]))*((1.+zarray)/chi)*dndz*pk[xx]*(bz[0]+bz[1]*(1.+zarray)^2)
   c_l[i]=int_tabulated(zarray,fz,/double)
ENDFOR

;MAD Multiply in all the constants
mod_cl=(3./2.)*((h0*100.)^2.)*omega_m*(1./(c^2.))*c_l
mod_ell=lvals

IF keyword_set(outfile) THEN BEGIN
   openw,1,outfile
   FOR i=0L,n_elements(lvals)-1 DO BEGIN
      printf,1,lvals[i],mod_cl[i]
   ENDFOR
   close,1
ENDIF

IF keyword_set(plotfile) THEN BEGIN
   PS_start,filename=plotfile,xsize=8,ysize=6
   plot,lvals,mod_cl*(10.^6.),linestyle=1,xtit='l',ytit='C_l x 10^6',/xlog,/ylog,xra=[6,2500],xsty=1,yra=[0.008,2],ysty=1
   PS_end,/png
ENDIF

return
END



