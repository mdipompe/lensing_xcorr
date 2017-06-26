;+
;  NAME:
;    xcorr_noise_spice
;  PURPOSE:
;    Cross correlate delta map with many different simulated noise
;    maps or rotated maps, using PolSpice (xcorr_noise uses anafast)
;
;  USE:
;    xcorr_noise_spice,delta,mask,bins=binsize,outfile=outfile,noise_loc='directory'
;
;  INPUT:
;    delta     - a HEALPix map of delta values
;    mask      - a HEALPix map of your mask (0 bad pixel, 1 good)
;    outfile   - string name of the output for all iterations
;                Default to 'xcorr_noise_binned.txt'
;    noise_loc - string name of directory where noise maps live.
;                Default to current directory. Looks for all fits
;                files in this directory.
;    bins      - binning in ell. Should match that used for the
;                measurement. Defaults to 5 per dex.
;
;  KEYWORDS:
;    maxell - maximum ell value (defaults to 3000)
;    
;  OUTPUT:
;    Placed in output file
;
;  HISTORY:
;    6-24-17 - Written - MAD (Dartmouth)
;-
PRO xcorr_noise_spice,delta,mask,noise_loc=noise_loc,bins=bins,outfile=outfile,maxell=maxell

IF ~keyword_set(noise_loc) THEN noise_loc='./'
IF ~keyword_set(outfile) THEN outfile='xcorr_noise_binned.txt'
IF ~keyword_set(maxell) THEN maxell=3000
IF ~keyword_set(bins) THEN bins=5

;MAD If output file already exists, don't run
check=file_search(outfile)
IF (check NE '') THEN BEGIN
   print,'XCORR_NOISE: Output file already exists, returning...'
   return
ENDIF

;MAD Get list of noise map fits files
rot_maps=file_search(noise_loc+'*.fits')
IF (rot_maps[0] EQ '') THEN message,'No maps found, exiting.'

openw,lun,outfile,/get_lun
printf,lun,';bin edge, bin_cent, bin width, c_l'

IF (n_elements(bins) EQ 1) THEN $
   bin_edges=define_ell_bins(bins,maxell) ELSE $
      bin_edges=bins
widths=bin_edges[1:n_elements(bin_edges)-1]-bin_edges[0:n_elements(bin_edges)-2]

;MAD Begin loop over noise maps
FOR i=0L,n_elements(rot_maps)-1 DO BEGIN
   counter,i,n_elements(rot_maps)

   spice4idl,mapfile1=delta,mapfile2=rot_maps[i],maskfile1=mask,$
             maskfile2=mask,clfile='raw_cl_tmp.txt',$
             bins=bins,binnedout='cl_binned_tmp.txt',$
             nlmax=maxell
   
   ;MAD Read binned cls, write out to single file
   readcol,'cl_binned_tmp.txt',binnedell,binnedcl,format='D'
   FOR j=0L,n_elements(binnedell)-1 DO BEGIN
      printf,lun,bin_edges[j],binnedell[j],widths[j],binnedcl[j],i,format='(F,1x,F,1x,F,1x,E,1x,I)'
   ENDFOR
   cmd=['rm','raw_cl_tmp.txt','cl_binned_tmp.txt']
   spawn,cmd,/noshell
ENDFOR
close,lun

return
END
