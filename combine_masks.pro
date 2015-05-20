;+
;  NAME:
;    combine_masks
;  PURPOSE:
;    Combine two binary masks
;
;  USE:
;    combine_masks,mask1,mask2,outfile=outfile
;
;  INPUT:
;    mask1 - a HEALPix map of binary mask
;    mask2 - a HEALPix map of binary mask 2
;
;  Optional Inputs:
;    outfile - string name of fits output
;
;  OUTPUT:
;    mask - HEALPix map of combined bitmask
;
;  HISTORY:
;    11-7-14 - Written - MAD (UWyo)
;-
PRO combine_masks,mask1,mask2,mask,outfile=outfile

;MAD If output file already exists, don't run just read it in
IF ~keyword_set(outfile) THEN check='' ELSE check=file_search(outfile)
IF (check NE '') THEN BEGIN
   print,'COMBINE_MASKS: Output file already exists, reading in and returning...'
   read_fits_map,outfile,mask
   return
ENDIF

;MAD Combine masks, keeping only pixels used in both
mask=intarr(n_elements(mask1))
mask[where((mask1 EQ 1) AND (mask2 EQ 1))]=1

;MAD Write out if needed
IF keyword_set(outfile) THEN write_fits_map,outfile,mask,/nest

return
END
