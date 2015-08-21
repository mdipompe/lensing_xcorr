;+
;  NAME:
;    combine_camb
;  PURPOSE:
;    Takes output files from CAMB and puts them together into a
;    power spectrum as a function of z, k, and ell
;
;  USE:
;   combine_camb,'path',zin,outstruct,maxell=maxell,chiin=chiin,outfile=outfile
;
;  INPUT:
;   path - string path to location of output files from
;          CAMB. Those files should have "matterpower" in their file
;          name. (Give trailing /)
;   zin - list of z values CAMB calculated power spectrum at.  
;
;  OPTIONAL INPUT
;   maxell - if set, will reinterpolate to a uniform grid in ell
;            (l=1,2,3,4...) up to maxell.  
;   chiin - if maxell is set, supply comoving distances at each z (in
;           Mpc/h) so k can be converted to ell.
;
;  OUTPUT:
;   outstruct - name of a structure for output with tags k, z, pk, and
;              l if maxell is set.
;   outfile - if set, writes structure out to file
;
;  HISTORY:
;    11-16-14 - Written - MAD (UWyo)
;     8-12-15 - No longer automatically interpolates to even grid in
;               ell, unless maxell is set.  If not, will just keep the
;               k spacing that camb used, and combine all files into
;               single fits file - MAD (UWyo)
;-
PRO combine_camb,path,zin,outstruct,maxell=maxell,chiin=chiin,outfile=outfile

;MAD If output file already exists, don't run just read it in
check=file_search(outfile)
IF (check NE '') THEN BEGIN
   print,'COMBINE_CAMB: Output file already exists, reading in and returning...'
   outstruct=mrdfits(outfile,1)
   return
ENDIF

;MAD Find CAMB output files
files=file_search(path+'*matterpower*')

;MAD Make sure there are as many z values as CAMB files.
;Doesn't guarantee that they match correctly, but better than no check
IF (n_elements(files) NE n_elements(zin)) THEN BEGIN
   print,'COMBINE_CAMB: Different number of z values than CAMB files, I quit.'
   return
ENDIF

;MAD interpolate to ell, if needed
IF keyword_set(maxell) THEN BEGIN
   lvals=findgen(maxell)+1
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
ENDIF ELSE BEGIN
   FOR i=0L,n_elements(files)-1 DO BEGIN
      readcol,files[i],kvals,pktmp,format='D',/silent
      tempz=fltarr(n_elements(kvals))+zin[i]
      IF (n_elements(k) EQ 0) THEN k=kvals ELSE k=[k,kvals]
      IF (n_elements(pk) EQ 0) THEN pk=pktmp ELSE pk=[pk,pktmp]
      IF (n_elements(zout) EQ 0) THEN zout=tempz ELSE zout=[zout,tempz]
   ENDFOR
   outstruct={k:0., z:0., pk:0.}
   outstruct=replicate(outstruct,n_elements(k))
   outstruct.k=k
   outstruct.z=zout
   outstruct.pk=pk
ENDELSE

IF keyword_set(outfile) THEN mwrfits,outstruct,outfile,/create

return
END

