;+
;  NAME:
;    rotate_map
;  PURPOSE:
;    Rotate a HEALPix map in galactic longitude (can also do
;    reflection in latitude)
;    
;  USE:
;   rotate_map,map,ang,rotmap,outfile='outfile.fits',/lat,/plot
;
;  INPUT:
;   map - HEALPix map to rotate
;   ang - rotation amount in longitude (degrees)
;
;  OPTIONAL INPUT:
;   outfile - string name of output file
;
;  KEYWORDS:
;   plot - if set, plots mollweide projections of original and rotated
;          maps (in galactic coordinates)
;   lat - if set, also rotate 180 degrees in latitude
;
;  OUTPUT:
;   rotmap - the rotated map
;
;  HISTORY:
;    11-21-14 - Written - MAD (UWyo)
;-
PRO rotate_map,map,ang,rotmap,outfile=outfile,plot=plot,lat=lat

;MAD If output file already exists, don't run just read it in
IF ~keyword_set(outfile) THEN check='' ELSE check=file_search(outfile)
IF (check NE '') THEN BEGIN
   print,'ROTATE_MAP: Output file already exists, reading in and returning...'
   read_fits_map,outfile,rotmap
   return
ENDIF

;MAD Plot original map if needed
IF keyword_set(plot) THEN $
   mollview,map,title='Original',/nest,max=max(map),min=min(map),coord='G',window=0

;MAD Determine nside, generate pixel numbers, coordinates
nside=SQRT(n_elements(map)/12.)
pix=indgen(12.*nside^2.,/long)
pix2ang_nest,nside,pix,theta,phi
galra=phi*(180./!dpi)
galdec=90.-((180./!dpi)*theta)
IF keyword_set(lat) THEN theta=((galdec*(-1))-90.)*(!dpi/((-1.)*180.))

;MAD Rotate
newgalra=galra+(ang)
xx=where(newgalra GE 360,cnt)
IF (cnt NE 0) THEN newgalra[xx]=newgalra[xx]-360.
newphi=newgalra*(!dpi/180.)
ang2pix_nest,nside,theta,newphi,newpix
match2,pix,newpix,pixm,newpixm
rotmap=map[newpixm]

IF keyword_set(plot) THEN BEGIN
   IF ~keyword_set(lat) THEN $
      titstring='Rotated Map, '+strtrim(ang,2)+' deg long, 0 deg lat' ELSE $
         titstring='Rotated Map, '+strtrim(ang,2)+' deg long, 180 deg lat'
   mollview,rotmap,title=titstring,/nest,max=max(map),min=min(map),coord='G',window=1
ENDIF

IF keyword_set(outfile) THEN write_fits_map,outfile,rotmap,/nest

return
END
