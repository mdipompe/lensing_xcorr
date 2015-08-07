;+
;  NAME:
;    alm_2_map
;  PURPOSE:
;    convert a healpix alm map to a healpix map (in nested format)
;
;  USE:
;    alm_2_map,'almfile',outmap,nside=nside,lmax=lmax
;
;  INPUT:
;    almfile - string name of file with alms
;
;  OPTIONAL INPUT:
;    nside - set nside of output map.  Defaults to 2048
;    lmax - set maximum multipole.  Defaults to 2048
;    outfile - set to output map to fits file
;
;  KEYWORDS:
;    full - if set, only includes pixels that have full area
;    
;  OUTPUT:
;    mask - HEALPix map of binary mask
;
;  HISTORY:
;    5-1-2015 - Written - MAD (UWyo)
;-
PRO alm_2_map,almfile,outmap,nside=nside,lmax=lmax,outfile=outfile

IF ~keyword_set(nside) THEN nside=2048
IF ~keyword_set(lmax) THEN lmax=2048

isynfast,'',outmap,alm_in=almfile,nside=nside,lmax=lmax
outmap=reorder(outmap[*,0],/r2n)

IF keyword_set(outfile) THEN $
   write_fits_map,outfile,outmap,/nest

END
