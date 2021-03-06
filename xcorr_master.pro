;+
;  NAME:
;    xcorr_master
;  PURPOSE:
;    Wrapper for all cross-correlation codes.  
;
;  USE:
;   xcorr_master,'planck_file','area_file','data_file'
;
;  INPUT:
;    planck_file - string name of file from planck containing kappa
;                  map.
;    planck_mask - string name of file with planck mask.
;    area_file - HEALPix map of pixel areas (from pix_area.pro)
;    data_file - file containing data, with ra and dec tags (in
;                equatorial coordinates)
;    sim_loc - path to simulated or rotated maps (with trailing /)
;
;  OPTIONAL INPUT:
;    binning - integer value of bins per dex (3, 4, or 5).  Defaults
;              to 5.    
;    ras - if plots is set, needs this to know plot range.  2 element
;          array of ra range
;    decs - if plots is set, needs this to know plot range.  2 element
;           array of dec range
;    h0 - little h.  Defaults to 0.702
;    omega_m - Omega_matter.  Defaults to 0.273
;    omega_b - Omega baryon.  Defaults to 0.046
;    omega_l - Omega_lambda.  Defaults to 0.727
;    min_ell - minimum ell to fit the bias (Defaults to 10)
;    max_ell - maximum ell to fit the bias (Defaults to 1000)
;
;  KEYWORDS:  
;    usepartial - use pixels that overlap mask components.  Default is
;                 to only use complete pixels
;    plots - if set, make nice plots of overlaid and stacked maps
;
;  OUTPUT:
;    
;  HISTORY:
;    11-16-14 - Written - MAD (UWyo)
;     7-16-15 - Modified to use CAMB4IDL - MAD (UWyo)
;     8-21-15 - Simplified model generation - MAD (UWyo)
;-
PRO xcorr_master,planck_file,planck_mask,areafile,datafile,sim_loc,usepartial=usepartial,$
                 binning=binning,plots=plots,ras=ras,decs=decs,h0=h0,omega_m=omega_m,$
                 omega_b=omega_b,omega_l=omega_l,min_ell=min_ell,max_ell=max_ell

;MAD Set start time
st=systime(1)

;MAD Set default binning
IF ~keyword_set(binning) THEN binning=5

;MAD Set default cosmology
IF ~keyword_set(h0) THEN h0=0.702
IF ~keyword_set(omega_m) THEN omega_m=0.275
IF ~keyword_set(omega_l) THEN omega_l=0.725

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;MAD Start calculating cross-power
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;MAD Read in planck data
read_fits_map,planck_file,kappa
read_fits_map,planck_mask,planckmask

;MAD Read in area file, build mask
read_fits_map,areafile,mapareas
IF keyword_set(usepartial) THEN $
   area_to_mask,mapareas,datamask,outfile='data_mask.fits' ELSE $
      area_to_mask,mapareas,datamask,/full,outfile='data_mask.fits'

;MAD Combine lensing and data masks into one
combine_masks,planckmask,datamask,mask,outfile='mask.fits'

;MAD Get density of data in pixels, make delta map
data=mrdfits(datafile,1)
pix_density,data,mapareas,mask,delta,outroot='data'

;MAD Cross correlate the maps
cross_corr,kappa,delta,mask,cl,ell,bins=binning,outroot='cl'


;MAD Rotate planck maps if needed
angles=20
WHILE (max(angles) LT 340) DO BEGIN
   angles=[angles,max(angles)+20]
ENDWHILE

FOR i=0L,n_elements(angles)-1 DO BEGIN
   outname='rotated/planck_rot_'+strtrim(angles[i],2)+'.fits'
   rotate_map,kappa,angles[i],outmap,outfile=outname
ENDFOR
FOR i=0L,n_elements(angles)-1 DO BEGIN
   outname='rotated/planck_rot2_'+strtrim(angles[i],2)+'.fits'
   rotate_map,kappa,angles[i],outmap,outfile=outname,/lat
ENDFOR


;MAD Cross correlate with simulated/rotated maps
xcorr_noise,delta,mask,noise_loc=sim_loc,bins=binning,/raw,root='sim_cl'

;MAD Calculate errors from cross-corrs with simulated/rotated maps
cross_corr_errs,'sim_cl_binned.txt',errs,covarmatrix,root='cl'

;MAD If wanted, make nice plots
IF keyword_set(plots) THEN BEGIN
   ;MAD Smooth delta map, reduce resolution
   gauss_smooth_map,delta,mask,60.,sdelta,outfile='smoothed_delta.fits'
   change_res,1024,sdelta,mask,sdelta1024,mask1024,outroot='delta'

   ;MAD Smooth kappa map
   gauss_smooth_map,kappa,mask,60.,skappa,outfile='smoothed_kappa.fits'
   change_res,1024,skappa,mask,skappa1024,mask1024,outroot='kappa'

   ;MAD Plot contour maps on each other
   kappa_delta_map,skappa1024,sdelta1024,mask1024,ras[0],ras[1],decs[0],decs[1],'kappa_delta_map.eps'

   ;MAD Bin and stack maps, make plot
   stack,sdelta,skappa,mask,deltstack,kapstack,outfile='stacked_data.eps'

   ;MAD Free up memory
   undefine,deltstack
   undefine,kapstack
   undefine,sdelta
   undefine,skappa
   undefine,sdelta1024
   undefine,skappa1024
   undefine,mask1024
ENDIF

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;MAD Generate model cross power
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
model_cross_corr,mod_ell,mod_cl,omega_m=omega_m,omega_l=omega_l,h0=h0,omega_b=omega_b,$
                 outfile='model_power.txt',plotfile='model_power.png'


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;MAD Fit cross-power, get bias
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
IF ~keyword_set(min_ell) THEN min_ell=10.
IF ~keyword_set(max_ell) THEN max_ell=1000.
fit_crosspower,cl,ell,min_ell,max_ell,mod_ell,mod_cl,bias,biaserr,err=errs,covar=covarmatrix,plotout='bias_fit.png'

;MAD Get total time
et=systime(1)
elapsed=et-st
print,'XCORR_MASTER: Total time = '+strtrim(elapsed/60.,2)+' minutes'

return
END
