;+
;  NAME:
;    gauss_smooth_map
;  PURPOSE:
;    Apply Gaussian smoothing to a HEALPix map.  Assumed to be in
;    nested format
;
;  USE:
;    gauss_smooth_map,map,mask,fwhm_arcmin,mapout,outfile='smoothedmap.fits',missingvals=missingvals,/plotmap
;
;  INPUT:
;    map - input HEALPix map to smooth
;    mask - binary mask to apply to data
;    fwhm - FWHM of Gaussian to apply in arcmin
;
;  OPTIONAL INPUT:
;    outfile - string name of file to write smoothed map to
;    missingvals - value to give to masked pixels (defaults to 0)
;
;  OPTIONAL KEYWORDS:
;    plotmap - if set will use mollview to plot the map
;
;  OUTPUT:
;    mapout - the smoothed HEALPix map
;
;  HISTORY:
;    11-11-14 - Written - MAD (UWyo)
;-
PRO gauss_smooth_map,map,mask,fwhm,outmap,missingvals=missingvals,outfile=outfile,plotmap=plotmap

;MAD If output file already exists, don't run just read it in
IF ~keyword_set(outfile) THEN check='' ELSE check=file_search(outfile)
IF (check NE '') THEN BEGIN
   print,'GAUSS_SMOOTH_MAP: Output file already exists, reading in and returning...'
   read_fits_map,outfile,outmap
   return
ENDIF

;MAD Set missingval to 0 if not supplied, mask data
IF ~keyword_set(missingval) THEN missingval=0
map[where(mask EQ 0)]=missingval

;MAD Smooth the map
ismoothing,map,outmap,fwhm_arcmin=fwhm,/nest

;MAD Plot if needed
IF keyword_set(plotmap) THEN $
   mollview,outmap,title='smoothed map',/nest,max=max(outmap),min=min(outmap)

;MAD write out smoothed map if desired
IF keyword_set(outfile) THEN write_fits_map,outfile,outmap,/nest


return
END
