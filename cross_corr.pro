;+
;  NAME:
;    cross_corr
;  PURPOSE:
;    cross-correlate two healpix maps
;
;  USE:
;    cross_corr,map1,map2,mask,cl,ell,bins=n_bins_per_dex,
;               outroot='out_file_root'
;
;  INPUT:
;    map1 - HEALPix map
;    map2 - HEALPix map
;    mask - binary HEALPix mask
;
;  OPTIONAL INPUT:
;    outroot - string name of the root of the file names (will be 
;               appended with _raw.txt, etc).  Defaults to 'cl'
;    bins_per_dex - integer number of bins per dex (3,4, or 5).
;                   Defualts to 5.
;    maxl - the maximum ell value to use in harmonic transform calculation
;
;    
;  OUTPUT:
;    cl - the cross-correlation coefficients
;    ell - the l values of the bin centers
;
;  HISTORY:
;    11-8-14 - Written - MAD (UWyo)
;     6-9-17 - Added maxl option - MAD (Dartmouth)
;-
PRO cross_corr,map1,map2,mask,cl,ell,bins=bins,outroot=outroot,maxl=maxl

IF ~keyword_set(bins) THEN bins=5
IF ~keyword_set(maxl) THEN maxl=3000

;MAD If output file already exists, don't run just read it in
IF ~keyword_set(outroot) THEN check='' ELSE check=file_search(outroot+'_binned_'+strtrim(bins,2)+'.txt')
IF (check NE '') THEN BEGIN
   print,'CROSS_CORR: Output file already exists, reading in and returning...'
   readcol,outroot+'_binned_'+strtrim(bins,2)+'.txt',edge,ell,width,cl,format='F'
   return
ENDIF

;MAD Shift means of maps to 0
IF (mean(map1[where(mask EQ 1)]) GT $
    stddev(map1[where(mask EQ 1)])) THEN BEGIN
   read,'XCORR_NOISE: Mean of map 1 more than 1-sigma from 0...continue? (y=0/n=1) ',input
   IF (input EQ 1) THEN message,'Problem with map 1, exiting.'
ENDIF
map1=map1-mean(map1[where(mask EQ 1)])

IF (mean(map2[where(mask EQ 1)]) GT $
    stddev(map2[where(mask EQ 1)])) THEN BEGIN
   read,'XCORR_NOISE: Mean of map 2 more than 1-sigma from 0...continue? (y/n) ',input
   IF (strtrim(input,2) NE 'y') THEN message,'Problem with map 2, exiting.'
ENDIF
map2=map2-mean(map2[where(mask EQ 1)])


;MAD deterime area correction factor
factor=n_elements(where(mask EQ 1))*(1./n_elements(mask))


;MAD Run ianafast to get cross-correlation
ianafast,map1,crosscorr,map2_in=map2,maskfile=mask,/nested,nlmax=maxl
crosscorr=crosscorr*(1./factor)

;MAD If needed, write out raw c_l values
IF keyword_set(outroot) THEN BEGIN
   openw,1,outroot+'_raw.txt'
   FOR i=0L,n_elements(crosscorr)-1 DO BEGIN
      printf,1,crosscorr[i]
   ENDFOR
   close,1
ENDIF

;MAD Define edges of bins (starting at l=10)
IF (bins EQ 3) THEN step=0.333
IF (bins EQ 4) THEN step=0.25
IF (bins EQ 5) THEN step=0.2
power=1.
bin_edge=0
WHILE (max(bin_edge) LT maxl) DO BEGIN
  IF (power EQ 1.) THEN bin_edge=10.^power ELSE bin_edge=[bin_edge,10.^power]
  power=power+step
ENDWHILE


;MAD Bin in l space
bin_llcl,crosscorr,bin_edge,ell,cl,deltal=bin_width,/uniform

;MAD If needed, write out binned data
IF keyword_set(outroot) THEN BEGIN
   openw,1,outroot+'_binned_'+strtrim(bins,2)+'.txt'
   printf,1,';bin_edge, bin_cent,bin_width,c_l'
   FOR i=0L,n_elements(cl)-1 DO BEGIN
      printf,1,bin_edge[i],ell[i],bin_width[i],cl[i],format='(F,1x,F,1x,F,1x,E)'
   ENDFOR
   close,1
ENDIF

return
END
