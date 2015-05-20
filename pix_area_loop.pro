;+
;  NAME:
;    pix_area_loop
;  PURPOSE:
;    Calculate area of HEALPix pixels in presence of mask with holes,
;    with control over the density per pixel desired for the area
;    calculation.
;
;    Relies on using Mangle representations of larger-area (compared
;    to your HEALPix maps, up to nside=64) HEALPix pixels.  Ransacks
;    the larger pixel, applies mask, and then calculates remaining
;    area in smaller sub-pixels.  More efficient than trying to ransack
;    the full area you use, especially if it is large and you need
;    high density. Works because of nested nature of HEALPix
;
;    Make sure your mask is in the same coordinate system as the HEALPix
;    pixels! (Typically Galactic)
;
;  USE:
;    pix_area_loop,manglepix,mask,nside,area_map,N=N,outfile='area_map.fits',scheme='6s'
;
;  INPUT:
;    manglepix - the Mangle HEALPix polygon file
;    mask - Mangle polygon mask describing holes in data
;    nside - nside of desired pixel area map for output
;
;  OPTIONAL INPUT:
;    N - number of points needed in each small pixel for area
;        calculation.  Defaults to 100.
;    outfile - filename for output, (fits map)
;    scheme - pixelization scheme of the mask (defaults to 6s)
;
;  OUTPUT:
;    areamap - HEALPix map containing area of each pixel
;
;  NOTES:
;    requires ransack.pro, an IDL wrapper for the Mangle tool ransack
;
;  HISTORY:
;    5-20-15 - Written - MAD (UWyo)
;-
PRO pix_area_loop,pixels,mask,nside,area_map,N=N,outfile=outfile,scheme=scheme

;MAD Set pixelization scheme
IF ~keyword_set(scheme) THEN scheme='6s'

;MAD Determine nside of larger pixels
nside_large=SQRT(n_elements(pixels)/12.)

;MAD This is a hack because the mangle representations of HEALPix 
;MAD pixels aren't perfect. In some cases ransacked points will lie 
;MAD slightly outside the actual pixel boundry.  When looping later, 
;MAD will only keep areas for high resolution pixels with centers 
;MAD actually inside the low resolution one.
;MAD Generate pixel numbers of higher resolution, find theta/phi
nums=lindgen(12.D*(nside^2.))
pix2ang_nest,nside,nums,theta,phi
;Find which larger pixel each of those smaller pixels falls in (nums)
ang2pix_nest,nside_large,theta,phi,nums

;MAD Set full pixel areas
pixarea=((4.*!dpi)*(180./!dpi)^2.)*(1./n_elements(nums))
pixarea_large=((4.*!dpi)*(180./!dpi)^2.)*(1./n_elements(pixels))

;MAD Calculate number of points needed in larger pixels to reach
;MAD density in smaller pixels
IF ~keyword_set(N) THEN N=100.
N_large=ceil((N/pixarea)*pixarea_large)

;MAD Loop over low resolution pixels
FOR i=0L,n_elements(pixels)-1 DO BEGIN
   st=systime(1)
   print,'Working on '+strtrim(i+1,2)+' of ',strtrim(n_elements(pixels),2)

   ;MAD Limit mask components to those overlapping pixel (for speed)
   where_polygons_overlap,pixels[i],mask,m1,nmatches
   IF (m1[0] NE -1) THEN BEGIN
      usemask=mask[m1]
      ;MAD Write temp polygon for ransack to use
      write_mangle_polygons,'tmp.ply',pixels[i]
      ransack,N_large,'tmp.ply','tmp.rsack',seed=616

      ;MAD Apply mask to ransacked data, calculate area used
      apply_poly_mask,'tmp.rsack',usemask,rin,rout,scheme=scheme,coords='Q'
      euler,rin.ra,rin.dec,ra1,dec1,2
      rin.ra=ra1
      rin.dec=dec1
      euler,rout.ra,rout.dec,ra1,dec1,2
      rout.ra=ra1
      rout.dec=dec1
      pix_area,2048,rin,rout,area_tmp
   ENDIF ELSE BEGIN
      area_tmp=dblarr(n_elements(nums))+pixarea
   ENDELSE

   ;MAD Use hack above - set area of any small pixel outside
   ;MAD large pixel boundry to 0.
      xx=where(nums NE i)
      area_tmp[xx]=0.

   ;MAD Add areas to final output area map
   IF (n_elements(area_map) EQ 0) THEN BEGIN
      area_map=area_tmp
   ENDIF ELSE BEGIN
      xx=where(area_tmp NE 0,cnt)
      IF (cnt NE 0) THEN BEGIN
         yy=where(area_map[xx] NE 0,cnt)
         IF (cnt NE 0) THEN BEGIN
            print,'Pixel area calculated twice, something''s wrong!'   ;Hack above failed,try again.
            stop
         ENDIF
         area_map=area_map+area_tmp
      ENDIF
   ENDELSE
et=systime(1)
print,'Loop took '+strtrim((et-st)*60.,2)+' minutes to run.'
ENDFOR

;MAD Write area map in healpix format, if wanted
IF keyword_set(outfile) THEN write_fits_map,outfile,area_map,/nest

return
END
