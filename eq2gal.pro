;+
;  NAME:
;    eq2gal
;  PURPOSE:
;    Take a healpix map in equatorial coordinates and convert it into
;    galactic. Assumes nested structure.
;
;  USE:
;   eq2gal,map,newmap,outfile='out.fits',/plot
;
;  INPUT:
;    map - HEALPix map to change coordinates of
;
;  OPTIONAL KEYWORDS: 
;    plot - set to make mollweide projection plots of in and out maps
;
;  OPTIONAL INPUT
;    outfile - string name of output fits map, if desired
;    
;  OUTPUT:
;    newmap - the new map
;
;  HISTORY:
;    2-20-15 - Written - MAD (UWyo)
;-
PRO eq2gal,map,newmap,outfile=outfile,plot=plot

IF keyword_set(plot) THEN mollview,map,/nest,title='Original Map'

;MAD Determine nside, get pixel coordinates
nside=SQRT(n_elements(map)/12.)
pix=indgen(12.*nside^2.,/long)
pix2ang_nest,nside,pix,theta,phi
galra=phi*(180./!dpi)
galdec=90.-((180./!dpi)*theta)

;MAD Convert pixel positions
euler,galra,galdec,ra,dec,2
newphi=ra*(!dpi/180.)
newtheta=(90.-dec)*(!dpi/180.)

;Get new pixel numbers, assign map values
ang2pix_nest,nside,newtheta,newphi,newpix
match2,pix,newpix,pixm,newpixm
newmap=map[newpixm]

;MAD Plot new map (if plot set)
IF keyword_set(plot) THEN mollview,newmap,/nest,title='New map'

;MAD Write out new map to fits (if outfile set)
IF keyword_set(outfile) THEN $
   write_fits_map,outfile,newmap,/nest

END
