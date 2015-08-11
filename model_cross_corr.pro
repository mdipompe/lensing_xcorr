;+
;  NAME:
;    model_cross_corr
;  PURPOSE:
;    Generate a model DM cross-power spectrum (bias=1)
;
;  USE:
;    model_cross_corr,mod_ell,mod_cl,power_spec='power_spec.fits',dndz='dndz.txt',zarray='z_chi.txt',$
;                      omega_m=omega_m,omega_l=omega_l,h0=h0,$
;                      outfile='model_power.txt'
;
;  INPUT:
;
;  Optional Inputs:
;    power_spec - structure or string name of fits file 
;                 containing the matter power
;                 spectrum.  Needs tags pk (the power spectrum), k
;                 (wavenumber), l (angular wavenumber), z (redshift).
;                 Can be made from CAMB data using camb4idl and
;                 combine_camb.pro. Defaults to 'power_spec_camb.fits'
;    dndz - string name of text file with dndz.  Will fit a spline
;           function and write it out to dndz_fit.txt.  You should
;           check this!!  Defaults to dndz.txt
;    zsample - string name of file containing z and chi values, as
;              made by chi_list.pro.  Defaults to z_chi.txt
;    omega_m - Omega_matter, defaults to 0.273.  Be sure it is
;              consistent with what you used when you calculated the
;              power spectrum!
;    omega_l - Omega_lambda, defaults to 0.727.  See above...
;    h0 - little h (H0/100), defaults to 0.702.  See above...
;    outfile - if supplied, writes model power out to text file
;    plotfile - if supplied, makes a plot of the model
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
;-
PRO model_cross_corr,mod_ell,mod_cl,power_spec=power_spec,dndz=dndz,zsample=zsample,omega_m=omega_m,omega_l=omega_l,h0=h0,outfile=outfile,plotfile=plotfile

;MAD If output file already exists, don't run just read it in
IF ~keyword_set(outfile) THEN check='' ELSE check=file_search(outfile)
IF (check NE '') THEN BEGIN
   print,'MODEL_CROSS_CORR: Output file already exists, reading in and returning...'
   readcol,outfile,mod_ell,mod_cl,format='F,D'
   return
ENDIF

;MAD Set some defaults
IF ~keyword_set(power_spec) THEN power_spec='power_spec_camb.fits'
IF ~keyword_set(dndz) THEN dndz='dndz.txt'
IF ~keyword_set(zsample) THEN zsample='z_chi.txt'

IF ~keyword_set(h0) THEN h0=0.702
IF ~keyword_set(omega_m) THEN omega_m=0.273
IF ~keyword_set(omega_l) THEN omega_l=0.727

;MAD Read in power spectrum if given file name instead of structure
check=size(power_spec)
IF ((check[0] EQ 0) AND (check[1] EQ 7)) THEN $
   pspec=mrdfits(power_spec,1) ELSE $
      pspec=power_spec


;Put h factors in power spectrum
delsq=pspec.pk*(1./(2.*!dpi^2.))*pspec.k^3.
pkk=pspec.k*h0
pkl=pspec.l
pkz=pspec.z
pk=delsq*(2.*!dpi^2.)/(pkk^3.)
undefine,pspec

;MAD Read in comoving distances 
readcol,zsample,z,chi,format='D'

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
fit_dndz,zdist,z,dndz


;MAD Do the integral part of C_l
c_l=fltarr(n_elements(lvals))
FOR i=0L,n_elements(c_l)-1 DO Begin
   xx=where(pkl EQ lvals[i])
   fz=((d_cmb[0]-chi)/(d_cmb[0]))*((1.+z)/chi)*dndz*pk[xx]
   c_l[i]=int_tabulated(z,fz,/double)
ENDFOR

;MAD Multiply in all the constants
c=2.99792458e5
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



