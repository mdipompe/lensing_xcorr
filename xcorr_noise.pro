;+
;  NAME:
;    xcorr_noise
;  PURPOSE:
;    Cross correlate delta map with many different simulated noise
;    maps or rotated maps
;
;  USE:
;    xcorr_noise,delta,mask,bins=binsize,root=root,noise_loc='directory',/raw
;
;  INPUT:
;    delta     - a HEALPix map of delta values
;    mask      - a HEALPix map of your mask (0 bad pixel, 1 good)
;    root      - string name of the root of output files
;                Default to 'xcorr_noise'
;    noise_loc - string name of directory where noise maps live.
;                Default to current directory
;
;  KEYWORDS:
;    raw - If set, writes out raw unbinned xcorr file too
;    
;  OUTPUT:
;    root_name+'_binned'
;  HISTORY:
;    11-7-14 - Written - MAD (UWyo)
;-
PRO xcorr_noise,delta,mask,noise_loc=noise_loc,bins=bins,raw=raw,root=root

IF ~keyword_set(noise_loc) THEN noise_loc='./'
IF ~keyword_set(root) THEN root='xcorr_noise'

;MAD If output file already exists, don't run
IF ~keyword_set(root) THEN check='' ELSE check=file_search(root+'_binned.txt')
IF (check NE '') THEN BEGIN
   print,'XCORR_NOISE: Output file already exists, returning...'
   return
ENDIF

;MAD Get list of noise map fits files
rot_maps=file_search(noise_loc+'*.fits')
IF (rot_maps[0] EQ '') THEN message,'No maps found, exiting.'

;MAD Shift mean delta to 0 (error if shift too much)
IF (mean(delta[where(mask EQ 1)]) GT $
    stddev(delta[where(mask EQ 1)])) THEN BEGIN
   read,'XCORR_NOISE: Mean delta more than 1-sigma from 0...continue? (y/n) ',input
   IF (strtrim(input,2) NE 'y') THEN message,'Problem with delta map, exiting.'
ENDIF
delta=delta-mean(delta[where(mask EQ 1)])

;MAD Calculate area scaling factor
factor=n_elements(where(mask EQ 1))*(1./n_elements(mask))

;MAD Open some output files for writing if wanted
IF keyword_set(raw) THEN BEGIN
   openw,1,root+'_raw.txt'
ENDIF
openw,2,root+'_binned.txt'
printf,2,';bin_edge, bin_cent,bin_width,c_l'

;MAD Generate ell bins (3, 4 or 5 per dex)
IF ~keyword_set(bins) THEN bins=5
power=1.
bin_edge=0
IF (bins EQ 3) THEN step=0.333
IF (bins EQ 4) THEN step=0.25
IF (bins EQ 5) THEN step=0.2
WHILE (max(bin_edge) LT 3000.) DO BEGIN
   IF (power EQ 1.) THEN bin_edge=10.^power ELSE bin_edge=[bin_edge,10.^power]
   power=power+step
ENDWHILE


;MAD Begin loop over noise maps
FOR i=0L,n_elements(rot_maps)-1 DO BEGIN
   counter,i,n_elements(rot_maps)
   read_fits_map,rot_maps[i],kappa
   ;MAD Shift mean of kappa map to 0
   kappa=kappa-mean(kappa[where(mask EQ 1)])

   ;MAD Call HEALPix cross-corr code
   ianafast,delta,crosscorr,map2_in=kappa,maskfile=mask,/nested,nlmax=3000
   crosscorr=crosscorr*(1./factor)

   ;MAD Write out raw c_l values, if wanted
   IF keyword_set(raw) THEN BEGIN
      FOR j=0L,n_elements(crosscorr)-1 DO BEGIN
         printf,1,crosscorr[j],i
      ENDFOR
   ENDIF

   ;MAD Call HEALPix binning code, write out
   bin_llcl,crosscorr,bin_edge,bin_cent,binned_cc,deltal=bin_width,/uniform 
   FOR j=0L,n_elements(binned_cc)-1 DO BEGIN
      printf,2,bin_edge[j],bin_cent[j],bin_width[j],binned_cc[j],i,format='(F,1x,F,1x, F,1x,E,1x,I)'
   ENDFOR
ENDFOR
IF keyword_set(raw) THEN close,1
close,2




END
