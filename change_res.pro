;+
;  NAME:
;    change_res
;  PURPOSE:
;    Reduce the resolution of a HEALPix map.  Is always pessimistic
;    (reduced res pixels overlapping masked full res pixels are
;    assumed to be masked) 
;
;  USE:
;   change_res,newnside,map,mask,newmap,newmask,outroot='outroot'
;
;  INPUT:
;    newnside - the reduced nside
;    map - HEALPix map to reduce resolution of
;    mask - binary mask for map
;
;  OPTIONAL INPUT
;    missingval - the value for the map where data is masked.
;                 Defaults to 0.
;    outroot - string name of root of output files.  Will be appended
;              with _mask and _data, as well as new nside
;    
;  OUTPUT:
;    newmap - the new reduced resolution map
;    newmask - the new reduced resolution mask map 
;
;  HISTORY:
;    11-16-14 - Written - MAD (UWyo)
;-
PRO change_res,newnside,map,mask,newmap,newmask,outroot=outroot

;MAD If output file already exists, don't run just read it in
IF ~keyword_set(outroot) THEN check='' ELSE $
   check=file_search(outroot+'_mask_'+strtrim(round(newnside),2)+'.fits')
IF (check NE '') THEN BEGIN
   print,'CHANGE_RES: Output file already exists, reading in and returning...'
   read_fits_map,outroot+'_mask_'+strtrim(round(newnside),2)+'.fits',newmask
   read_fits_map,outroot+'_data_'+strtrim(round(newnside),2)+'.fits',newmap
   return
ENDIF

npixin=n_elements(map)
nsidein=SQRT(npixin/12.)

IF (nsidein EQ newnside) THEN BEGIN
   print,'CHANGE_RES: You aren''t changing the resolution! Skipping...'
   return
ENDIF

IF (newnside LT nsidein) THEN print,'CHANGE_RES: Reducing map resolution...'
IF (newnside GT nsidein) THEN print,'CHANGE_RES: Increasing map resolution...'

map[where(mask EQ 0)]=-99

ud_grade,mask,newmask,nside_out=newnside,order_in='nested',$
         order_out='nested',bad_data=0,/pessimistic

IF keyword_set(outroot) THEN BEGIN
   outfile=outroot+'_mask_'+strtrim(round(newnside),2)+'.fits'
   write_fits_map,outfile,newmask,/nested
ENDIF

ud_grade,map,newmap,nside_out=newnside,order_in='nested',$
         order_out='nested',bad_data=-99,/pessimistic

IF keyword_set(outroot) THEN BEGIN
   outfile=outroot+'_data_'+strtrim(round(newnside),2)+'.fits'
   write_fits_map,outfile,newmap,/nested
ENDIF

return
END
