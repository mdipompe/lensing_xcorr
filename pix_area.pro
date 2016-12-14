;+
;  NAME:
;    pix_area
;  PURPOSE:
;    Calculate area of HEALPix pixels in presence of mask with holes
;
;  USE:
;    pix_area,nside,inmask,outmask,areamap,outmap='outmap.fits'
;
;  INPUT:
;    nside - nside of the HEALPix map you want areas for
;    inmask - structure containing points inside mask holes (must
;             contain tags ra and dec)
;    outmask - structure containing points outside mask holes (must
;              contain tags ra and dec)
;             
;  OUTPUT:
;    areamap - HEALPix map containing area of each pixel
;
;  HISTORY:
;    11-16-14 - Written - MAD (UWyo)
;-
PRO pix_area,nside,inmask,outmask,areamap,outmap=outmap

;MAD If output file already exists, don't run just read it in
IF ~keyword_set(outmap) THEN check='' ELSE check=file_search(outmap)
IF (check NE '') THEN BEGIN
   print,'PIX_AREA: Output file already exists, reading in and returning...'
   read_fits_map,outmap,areamap
   return
ENDIF

IF (outmask[0].ra EQ -999. AND outmask[0].dec EQ -999.) THEN BEGIN
   outflag=1
   outmask[0].ra=0.
   outmask[0].dec=0.
ENDIF ELSE BEGIN
   outflag=0
ENDELSE
IF (inmask[0].ra EQ -999. AND inmask[0].dec EQ -999.) THEN BEGIN
   inflag=1
   inmask[0].ra=0.
   inmask[0].dec=0.
ENDIF ELSE BEGIN
   inflag=0
ENDELSE



;MAD Get area of each full healpix pixel
npix=12.*double(nside)^2.
pixarea=((4.*!dpi)*((180./!dpi)^2.))*(1./npix)

;MAD Convert RA and DEC to galactic, xyz
euler,outmask.ra,outmask.dec,galra_out,galdec_out,1
phi_out=galra_out*(!dpi/180.)
theta_out=(90.-galdec_out)*(!dpi/180.)
euler,inmask.ra,inmask.dec,galra_in,galdec_in,1
phi_in=galra_in*(!dpi/180.)
theta_in=(90.-galdec_in)*(!dpi/180.)

;MAD Get pixels of each position
ang2pix_nest,nside,theta_out,phi_out,ipix_out
ang2pix_nest,nside,theta_in,phi_in,ipix_in

;MAD Count number of objects in each pixel that are also inside and
;outside of mask components
n_out=double(histogram(ipix_out,binsize=1,min=(0),max=npix-1))
n_in=double(histogram(ipix_in,binsize=1,min=(0),max=npix-1))

IF (outflag EQ 1) THEN n_out=n_out*0.
IF (inflag EQ 1) THEN n_in=n_in*0.

;MAD Caclulate used area of each pixel
areamap=(n_out/(n_out+n_in))*pixarea
xx=where(n_out EQ 0,cnt)
IF (cnt NE 0) THEN areamap[xx]=0.

IF keyword_set(outmap) THEN write_fits_map,outmap,areamap,/nest

return
END
