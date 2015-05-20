;+
;  NAME:
;    convergence
;  PURPOSE:
;    Convert a Planck HEALPix lensing map to a lensing convergence map
;
;  USE:
;    convergence,'planck_lensing_map.fits',kappa,mask,$
;                 outroot='output_file_root'
;
;  INPUT:
;    file  - string name of file containing planck lensing map
;
;  OPTIONAL INPUT:
;     outroot - if set, outputs files with this name 
;               appended with _kappa.fits, etc).  
;
;    
;  OUTPUT:
;    kappa - HEALPix map of planck convergence
;    mask  - HEALPix map of planck mask
;
;  HISTORY:
;    11-8-14 - Written - MAD (UWyo)
;-
PRO convergence,file,kappa,mask,outroot=outroot

;MAD If output file already exists, don't run just read it in
IF ~keyword_set(outroot) THEN check='' ELSE check=file_search(outroot+'_kappa.fits')
IF (check NE '') THEN BEGIN
   print,'CONVERGENCE: Output file already exists, reading in and returning...'
   read_fits_map,outroot+'_kappa.fits',kappa
   read_fits_map,outroot+'_mask.fits',mask
   return
ENDIF


;MAD Read in the Plank potential map
read_fits_map,file, hplanck, hdr
mask=hplanck[*,1]

;MAD Read in the transfer function
tf_etc = mrdfits(file,2)

;MAD Convolve the map with the transfer function to get the
;lensing convergence map
p_ell = findgen(2049)
ismoothing, hplanck[*,0], kappa,  $
      beam=tf_etc[*].rlpp/(tf_etc[*].rlpp^2.+1.e6)*0.5*p_ell*(1.+p_ell), $
      /nest,lmax=3000

;MAD If wanted, write the maps to fits files
IF keyword_set(outroot) THEN BEGIN
      write_fits_map,outroot+'_kappa.fits',kappa,hdr,/nest
   write_fits_map,outroot+'_mask.fits',mask,hdr,/nest
ENDIF

return
END
