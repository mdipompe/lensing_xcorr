;+
;  NAME:
;    area_to_mask
;  PURPOSE:
;    Convert a map of pixel areas into a (binary)mask
;
;  USE:
;    area_to_mask,areas,/full,outfile='out.fits'
;
;  INPUT:
;    areas   - a HEALPix map of pixel areas
;
;  OPTIONAL INPUT:
;    outfile - string name of output file
;
;  KEYWORDS:
;    full - if set, only includes pixels that have full area
;    
;  OUTPUT:
;    mask - HEALPix map of binary mask
;
;  HISTORY:
;    11-7-14 - Written - MAD (UWyo)
;-
PRO area_to_mask,areas,mask,full=full,outfile=outfile

;MAD If output file already exists, don't run just read it in
IF ~keyword_set(outfile) THEN check='' ELSE check=file_search(outfile)
IF (check NE '') THEN BEGIN
   print,'AREA_TO_MASK: Output file already exists, reading in and returning...'
   read_fits_map,outfile,mask
   return
ENDIF

;MAD Keep any pixels that have some area used, or only those
;completely untouched if /full is set
mask=intarr(n_elements(areas))
IF ~keyword_set(full) THEN BEGIN
   mask[where(areas GT 0)]=1
ENDIF ELSE BEGIN
   mask[where(areas EQ max(areas))]=1
ENDELSE

;MAD Write out file, if needed
IF keyword_set(outfile) THEN write_fits_map,outfile,mask,/nest

return
END

