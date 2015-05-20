;+
;  NAME:
;    apply_poly_mask
;  PURPOSE:
;    Apply MANGLE polygon mask to random sample to figure out how much
;    area of each HEALPix pixel is used.
;
;  USE:
;    apply_poly_mask,'randfile',poly,randin,randout, $
;            outroot='root',coords='Q',scheme='6s'
;
;  INPUT:
;    randfile - string file of random positions (equatorial coords)
;    poly - MANGLE polygons to check in/out of.  Assumed in galactic coords
;
;  OPTIONAL INPUT:
;    outroot - if set, outputs fits files with "_out.fits" and
;              "_in.fits" appended of objects in/out of mask
;    coords - Coordinate system of polygons.  Defaults to galactic,
;             options are "Q" for equatorial, "G" for galactic, or "E"
;             for ecplitic
;    scheme - if polygons are pixelized, provide the scheme so
;             is_in_window_pix can speed things up
;
;  OUTPUT:
;    randin - structure of random points inside polygons
;    randout - structure of random points outside polygons
;
;  HISTORY:
;    11-7-14 - Written - MAD (UWyo)
;-
PRO apply_poly_mask,randfile,poly,randin,randout,outroot=outroot,coords=coords,scheme=scheme

;MAD Read in random objects
readcol,randfile,ra,dec

;MAD Convert random objects to system of polygons
IF keyword_set(coords) THEN BEGIN
   IF (coords EQ 'Q') THEN BEGIN
      long=ra
      lat=dec
   ENDIF ELSE IF (coords EQ 'G') THEN BEGIN
      euler,ra,dec,long,lat,1
   ENDIF ELSE IF (coords EQ 'E') THEN BEGIN
      euler,ra,dec,long,lat,3
   ENDIF
ENDIF ELSE BEGIN
   euler,ra,dec,long,lat,1
ENDELSE


;MAD Run is_in_window or is_in_window_pix
IF keyword_set(scheme) THEN $
   in=is_in_window_pix(ra=long,dec=lat,scheme=scheme,poly) ELSE $
      in=is_in_window(ra=long,dec=lat,poly)

;MAD Determine which ra and dec values are in/out of mask
xx=where(in EQ 1,cnt)
IF (cnt NE 0) THEN ra_inmask=ra[xx] ELSE ra_inmask=-999.
IF (cnt NE 0) THEN dec_inmask=dec[xx] ELSE dec_inmask=-999.
xx=where(in EQ 0,cnt)
IF (cnt NE 0) THEN ra=ra[xx] ELSE ra=-999.    ;Overwrites ra to save memory
IF (cnt NE 0) THEN dec=dec[xx] ELSE dec=-999. ;Overwrites dec to save memory

;MAD Set output structures
randout={ra:0.,dec:0.}
randout=replicate(randout,n_elements(dec))
randout.ra=ra
randout.dec=dec

randin={ra:0.,dec:0.}
randin=replicate(randin,n_elements(dec_inmask))
randin.ra=ra_inmask
randin.dec=dec_inmask

;MAD Write out files, if needed
IF keyword_set(outroot) THEN BEGIN
   mwrfits,randout,outroot+'_out.fits',/create
   mwrfits,randin,outroot+'_in.fits',/create
ENDIF

return
END
