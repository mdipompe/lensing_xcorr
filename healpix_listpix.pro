;+
;  NAME:
;    healpix_listpix
;  PURPOSE:
;    Take a healpix map and output list of pixel values and
;    coordinates (equatorial centers of pixels).  
;
;  USE:
;    healpix_listpix,map,ra,dec,rarange=[minra,maxra],decrange=[mindec,maxdec],outfile='out.txt'
;
;  INPUT:
;    map - input HEALPix map (Galactic coords, nested format)
;    
;  OPTIONAL INPUT:
;    outfile - string name of file to write pixel values and coords to
;    rarange - 2 element array with minimum and maximum RA values to
;              write out (decimal degrees)
;    decrange - same as above, for Dec.
;
;  OUTPUT:
;    ra, dec - coords of pixel centers (in same order as values in map)
;
;  HISTORY:
;    3-10-15 - Written - MAD (UWyo)
;-
PRO healpix_listpix,map,ra,dec,rarange=rarange,decrange=decrange,outfile=outfile

;MAD Find nside from healpix map
npix=n_elements(map)
nside=SQRT(npix/12.)

;MAD Generate positions of pixel centers
pixvals=findgen(npix)
pix2ang_nest,nside,pixvals,theta,phi
;MAD Convert to galactic
l=phi*(180./!dpi)
b=90.-(theta*(180./!dpi))
;MAD Convert to equatorial
euler,l,b,ra,dec,2

;Limit to RA and DEC ranges, if required
IF (keyword_set(rarange)) THEN BEGIN
   xx=where(ra GE rarange[0] AND ra LE rarange[1])
   map=map[xx]
   ra=ra[xx]
   dec=dec[xx]
ENDIF
IF (keyword_set(decrange)) THEN BEGIN
   xx=where(dec GE decrange[0] AND dec LE decrange[1])
   map=map[xx]
   ra=ra[xx]
   dec=dec[xx]
ENDIF

;MAD Write out text file, if required
IF (keyword_set(outfile)) THEN BEGIN
   openw,1,outfile
   printf,1,'ra       dec       value'
   FOR i=0L,n_elements(map)-1 DO BEGIN
      printf,1,strtrim(ra[i],2) + '    ' + strtrim(dec[i],2) + '    ' + strtrim(map[i],2),format='(A)'
   ENDFOR
   close,1
ENDIF

return
END
