;+
;  NAME:
;    combine_camb
;  PURPOSE:
;    Takes output files from CAMB and puts them together into a
;    power spectrum as a function of z, k, and ell
;
;  USE:
;   combine_camb,'path',zin,chiin,outstruct,maxell=maxell,outfile=outfile
;
;  INPUT:
;   path - string path to location of output files from
;          CAMB. Those files should have "matterpower" in their file
;          name. (Give trailing /)
;   zin - list of z values CAMB calculated power spectrum at.  
;   chiin - list of comiving distances corresponding to each z.  In
;           units of Mpc/h (so if chi is in Mpc make sure to multiply
;           by h!).
;
;  OPTIONAL INPUT
;   maxell - maximum value of ell.  Defaults to 3000.
;
;  OUTPUT:
;   outstruct - name of a structure for output with tags k, l, z, pk.
;   outfile - if set, writes structure out to file
;
;  HISTORY:
;    11-16-14 - Written - MAD (UWyo)
;-
PRO combine_camb,path,zin,chiin,outstruct,maxell=maxell,outfile=outfile

;MAD If output file already exists, don't run just read it in
check=file_search(outfile)
IF (check NE '') THEN BEGIN
   print,'COMBINE_CAMB: Output file already exists, reading in and returning...'
   outstruct=mrdfits(outfile,1)
   return
ENDIF

;MAD Find CAMB output files
files=file_search(path+'*matterpower*')

;MAD Set max ell
IF ~keyword_set(maxell) THEN lvals=findgen(3000)+1 ELSE $
   lvals=findgen(maxell-1.)+1

;MAD Make sure there are as many z values as CAMB files.
;Doesn't guarantee that they match correctly, but better than no check
IF (n_elements(files) NE n_elements(zin)) THEN BEGIN
   print,'COMBINE_CAMB: Different number of z values than CAMB files, I quit.'
   return
ENDIF

;MAD Loop over files, convert ell to k, interpolate, build pk

FOR i=0L,n_elements(files)-1 DO BEGIN
  readcol,files[i],tempk,tempp,format='D',/silent
  kvals=lvals*(1./chiin[i])
  newp=interpol(tempp,tempk,kvals)
  tempz=fltarr(n_elements(kvals))+zin[i]
  IF (n_elements(k) EQ 0) THEN k=kvals ELSE k=[k,kvals]
  IF (n_elements(ell) EQ 0) THEN ell=lvals ELSE ell=[ell,lvals]
  IF (n_elements(pk) EQ 0) THEN pk=newp ELSE pk=[pk,newp]
  IF (n_elements(zout) EQ 0) THEN zout=tempz ELSE zout=[zout,tempz]
ENDFOR

outstruct={k:0., l:0., z:0., pk:0.}
outstruct=replicate(outstruct,n_elements(k))
outstruct.k=k
outstruct.l=ell
outstruct.z=zout
outstruct.pk=pk

IF keyword_set(outfile) THEN mwrfits,outstruct,outfile,/create

return
END

