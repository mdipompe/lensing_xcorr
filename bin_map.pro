;+
;  NAME:
;    bin_map
;  PURPOSE:
;    Get the mean value of a HEALPix map in annuli around some set of points
;
;  USE:
;    bin_map,map,data,bincents,binned_vals,mask=mask,coords='G',bins=bins,region=region,outfile='out.txt',/loop
;
;  INPUT:
;    map - HEALPix map to bin
;    data - structure of data to bin around (tags ra/dec or l/b)
;
;  OPTIONAL INPUT:
;    bins - array of bin edges (in arcmin).  Defaults to 10' bins from 0 to 3 degrees
;    outfile - string name of text file to write bincents and
;              binned_vals to
;    mask - binary map of useable pixels in map
;    region - mangle polygon (single) defining region of interest.
;             Speeds things up by excluding pixels/data outside this
;             region.  In Equatorial coordinates!
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
;-
PRO bin_map,map,data,bincents,binned_vals,mask=mask,coords=coords,bins=bins,region=region,outfile=outfile,loop=loop

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

;MAD Limit pixels to region
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
IF (cnt NE 0) THEN remove,masked,pixnums,newmap,pix_l,pix_b

;MAD Set default bin edges
IF ~keyword_set(binedges) THEN binedges=findgen((6.*3.)+1)*10.
;MAD Generate bin centers
FOR i=0L,n_elements(binedges)-2 DO BEGIN
   val=binedges[i]+((binedges[i+1]-binedges[i])/2.)
   IF (n_elements(bincents) EQ 0) THEN bincents=val ELSE bincents=[bincents,val]
ENDFOR

;MAD Locate pixels within max bin size, find total map value in annuli
;while also counting number of pixels matched in each for mean calculation
k=0L
;MAD Set step size for loop, or just make one big step if /loop not set
IF keyword_set(loop) THEN step=10000. ELSE step=n_elements(pix_l)
;MAD n will be the number of pixels in each bin
n=dblarr(n_elements(bincents))
binned_tot=dblarr(n_elements(bincents))
WHILE (k LT n_elements(pix_l)) DO BEGIN
   IF keyword_set(loop) THEN counter,k,n_elements(pix_l)
   tempindx=lindgen(step)+(k)
   spherematch,data_l,data_b,pix_l[tempindx],pix_b[tempindx],max(binedges/60.),m1,m2,sep,maxmatch=0
   IF (m1[0] NE -1) THEN BEGIN
      FOR i=0L,n_elements(binedges)-2 DO BEGIN
         xx=where(sep GE binedges[i]/60. AND sep LT binedges[i+1]/60.,cnt)
         IF (cnt NE 0) THEN BEGIN
            n[i]=n[i]+n_elements(xx)
            IF (i EQ 0) THEN binned_step=total(newmap[tempindx[m2[xx]]]) ELSE $
               binned_step=[binned_step,total(newmap[tempindx[m2[xx]]])]
         ENDIF ELSE BEGIN
            IF (i EQ 0) THEN binned_step=0. ELSE $
               binned_step=[binned_step,0.]
         ENDELSE
      ENDFOR
      binned_tot=binned_tot+binned_step
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

;MAD Calculate mean from total and n
binned_vals=binned_tot*(1./n)

;MAD Write out file, if needed
IF keyword_set(outfile) THEN BEGIN
   openw,1,outfile
   FOR i=0L,n_elements(bincents)-1 DO BEGIN
      printf,1,bincents[i],binned_vals[i],format='(F,1x,E)'
   ENDFOR
   close,1
ENDIF

return
END
