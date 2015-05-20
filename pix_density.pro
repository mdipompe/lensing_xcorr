;+
;  NAME:
;    pix_density
;  PURPOSE:
;    Calculate the relative object density in HEALPix pixels
;
;  USE:
;    pix_density,data,areas,mask,delta,N=N,rho=rho,outroot='outroot'
;
;  INPUT:
;    data - Data structure with at least ra and dec tags
;    areas - HEALPix map containing used area of each pixel
;
;  OPTIONAL INPUT
;    outroot - if set, will write out maps for number, density, and
;              relative density of each pixels
;    
;  OUTPUT:
;    delta - HEALPix map of relative density in each pixel
;
;  OPTIONAL OUTPUT
;    N - HEALPix map of number of objects in each pixel
;    rho - HEALPix map of density in heach pixel
;
;  HISTORY:
;    11-16-14 - Written - MAD (UWyo)
;-
PRO pix_density,data,areas,mask,delta,N=N,rho=rho,outroot=outroot

;MAD If output file already exists, don't run just read it in
IF ~keyword_set(outroot) THEN check='' ELSE check=file_search(outroot+'_delta.fits')
IF (check NE '') THEN BEGIN
   print,'PIX_DENSITY: Output file already exists, reading in and returning...'
   read_fits_map,outroot+'_delta.fits',delta
   read_fits_map,outroot+'_N.fits',N
   read_fits_map,outroot+'_rho.fits',rho
   return
ENDIF

;MAD Set HEALPix params based input maps
npix=n_elements(areas)
nside=SQRT(npix/12.)

;MAD Get used pixels from mask
use_pix=where(mask EQ 1)

;MAD Convert to Galactic coords, get phi and theta on sphere
euler,data.ra,data.dec,galra,galdec,1
phi=galra*(!dpi/180.)
theta=(90-galdec)*(!dpi/180.)

;MAD Find which pixel each point is in
ang2pix_nest,nside,theta,phi,ipix

;MAD Count number in each pixel
N=histogram(ipix,binsize=1,min=(0),max=npix-1)

;MAD Calculate density averaged over the whole region
total_dense=total(N[use_pix])*(1./total(areas[use_pix]))

;MAD Calculate density in each pixel
rho=N*(1./areas)

;MAD Calculate delta
delta=(rho-total_dense)*(1./total_dense)

;MAD If needed, write files
IF keyword_set(outroot) THEN BEGIN
   write_fits_map,outroot+'_N.fits',N,/nest
   write_fits_map,outroot+'_rho.fits',rho,/nest
   write_fits_map,outroot+'_delta.fits',delta,/nest
ENDIF

return
END
