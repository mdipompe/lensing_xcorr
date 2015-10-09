;+
;  NAME:
;    bin_map
;  PURPOSE:
;    Get the mean value of a HEALPix map in annuli around some set of points
;
;  USE:
;    bin_map,map,data,bincents,binned_vals,mask=mask,coords='G',binedges=binedges,
;                region=region,outfile='out.txt',weights=weights/loop
;
;  INPUT:
;    map - HEALPix map to bin
;    data - structure of data to bin around (tags ra/dec or l/b)
;
;  OPTIONAL INPUT:
;    binedges - array of bin edges (in arcmin).  Defaults to 10' bins from 0 to 3 degrees
;    outfile - string name of text file to write bincents and
;              binned_vals to
;    mask - binary map of useable pixels in map
;    region - mangle polygon (single) defining region of interest.
;             Speeds things up by excluding pixels/data outside this
;             region.  In Equatorial coordinates!
;    weights - weight for each data point
;
;  KEYWORDS:        
;    coords - coordinates of data (E'Q'uatorial or 'G'alactic).
;            Defaults to 'Q'.  If 'Q' need tags ra and dec, if 'G',
;            need tags 'l' and 'b'.
;    loop - set to break up kappa map into chunks so that you
;           don't have memory issues trying to spherematch all
;           your data to all the pixels.  Useful for high-res HEALPix
;           maps and/or large data sets
;
;  OUTPUT:
;    bincents - the bin centers (degrees)
;    binned_vals - mean map values in annuli
;
;  NOTES:
;    Assumed that HEALPix map is in Galactic coordinates, and data are
;    in equatorial (will transform for you)
;
;
;  HISTORY:
;    6-15-15 - Written - MAD (UWyo)
;    7-10-15 - Improved speed, especially for large datasets - MAD
;              (UWyo)
;-
PRO bin_map,map,data,bincents,binned_vals,mask=mask,coords=coords,binedges=binedges,$
            region=region,outfile=outfile,loop=loop,weights=weights

st=timer()

;MAD Set default coordinates
IF ~keyword_set(coords) THEN coords='Q'

;MAD Transform data to Galactic if needed
IF (coords EQ 'Q') THEN BEGIN
   euler,data.ra,data.dec,data_l,data_b,1
ENDIF ELSE IF (coords EQ 'G') THEN BEGIN
   data_l=data.l
   data_b=data.b
ENDIF

;MAD Get npix, nside of map (round is to make it integer)
nside=round(SQRT(n_elements(map)/12.))
pixnums=lindgen(n_elements(map))

;MAD Get pixel center coordinates
healpix_coords,nside,pix_l,pix_b

;MAD Get pixel information of data
phi=data_l*(!dpi/180.)
theta=(90-data_b)*(!dpi/180.)
ang2pix_nest,nside,theta,phi,ipix
;Count number of data points in each pixel
;h=histogram(ipix,bin=1,min=0,max=max(pixnums))
   
;MAD Limit pixels to region of interest
IF keyword_set(region) THEN BEGIN
   euler,pix_l,pix_b,pix_ra,pix_dec,2
   in=is_in_window(ra=pix_ra,dec=pix_dec,region)
ENDIF

;MAD Mask pixels (or not)
IF (keyword_set(region) AND keyword_set(mask)) THEN BEGIN
   masked=where(mask EQ 0 OR in EQ 0, cnt)
ENDIF ELSE IF (keyword_set(region) AND ~keyword_set(mask)) THEN BEGIN
   masked=where(in EQ 0, cnt)
ENDIF ELSE IF (keyword_set(mask) AND ~keyword_set(region)) THEN BEGIN
   masked=where(mask EQ 0, cnt)
ENDIF ELSE IF (~keyword_set(mask) AND ~keyword_set(region)) THEN BEGIN
   cnt=0
ENDIF
newmap=map
IF (cnt NE 0) THEN remove,masked,pixnums,newmap,pix_l,pix_b;,h

;MAD Identify pixels with data points in them
;data_pix=where(h GT 0)

;MAD Set default bin edges
IF ~keyword_set(binedges) THEN binedges=findgen((6.*3.)+1)*10.
;MAD Generate bin centers
IF keyword_set(bincents) THEN undefine,bincents
FOR i=0L,n_elements(binedges)-2 DO BEGIN
   val=binedges[i]+((binedges[i+1]-binedges[i])/2.)
   IF (n_elements(bincents) EQ 0) THEN bincents=val ELSE bincents=[bincents,val]
ENDFOR


;MAD Initialize some arrays to fill
total_vals=dblarr(n_elements(bincents))
num=dblarr(n_elements(bincents))
binned_vals=dblarr(n_elements(bincents))
k=0L
;MAD Loop over chunks (if \loop set)
IF keyword_set(loop) THEN step=10000. ELSE step=n_elements(pix_l)
WHILE (k LT n_elements(pix_l)) DO BEGIN
   IF keyword_set(loop) THEN counter,k,n_elements(pix_l)
   tempindx=lindgen(step)+(k)
   spherematch,pix_l[ipix],pix_b[ipix],pix_l[tempindx],pix_b[tempindx],max(binedges/60.),m1,m2,sep,maxmatch=0
   IF (m1[0] NE -1) THEN BEGIN
      FOR i=0L,n_elements(binedges)-2 DO BEGIN
         xx=where(sep GE binedges[i]/60. AND sep LT binedges[i+1]/60.,cnt)
         IF (cnt NE 0) THEN BEGIN
            total_vals[i]=total_vals[i]+total(newmap[tempindx[m2[xx]]])
            num[i]=num[i]+n_elements(xx)
         ENDIF
      ENDFOR
   ENDIF
   k=k+step
   IF (k+10000 LT n_elements(pix_l)) THEN BEGIN
      step=10000
   ENDIF ELSE BEGIN
      IF (k+1000 LT n_elements(pix_l)) THEN BEGIN
         step=1000
      ENDIF ELSE BEGIN
         IF (k+100 LT n_elements(pix_l)) THEN BEGIN
            step=100
         ENDIF ELSE BEGIN
            IF (k+10 LT n_elements(pix_l)) THEN BEGIN
               step=10
            ENDIF ELSE BEGIN
               step=1
            ENDELSE
         ENDELSE
      ENDELSE
   ENDELSE
ENDWHILE

;MAD Find which annuli had pixels matched
use=where(total_vals NE 0)
;MAD Calculate mean per annuli
binned_vals[use]=total_vals[use]*(1./num[use])

;MAD Write out file, if needed
IF keyword_set(outfile) THEN BEGIN
   openw,1,outfile
   FOR i=0L,n_elements(bincents)-1 DO BEGIN
      printf,1,bincents[i],binned_vals[i],format='(F,1x,E)'
   ENDFOR
   close,1
ENDIF

et=timer(st=st,/fin,unit='h')

return
END
