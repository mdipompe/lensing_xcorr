;+
;  NAME:
;    cross_corr_errs.pro
;  PURPOSE:
;    Get errors on cross-corr using results from xcorr_noise.pro
;
;  USE:
;    cross_corr_errs,'noisefile.txt',errs,cmatrix,root='rootname'
;
;  INPUT:
;    infile - text file containing ell, cl, and run number of noise simulations.
;             format should be 1) bin edge, 2) bin center, 3) bin
;             width, 4) cl, 5) run number
;
;  OPTIONAL INPUT:
;    root - string name of the root of optional output files (appended
;           with _errs.txt and _covar.txt)
;
;  OUTPUT:
;    errs - the 1 sigma errors
;    covar - the covariance matrix
;
;  HISTORY:
;    11-11-14 - Written - MAD (UWyo)
;-
PRO cross_corr_errs,infile,errs,C,root=root

;MAD If output file already exists, don't run just read it in
IF ~keyword_set(root) THEN check='' ELSE check=file_search(root+'_errs.txt')
IF (check NE '') THEN BEGIN
   print,'CROSS_CORR_ERRS: Output file already exists, reading in and returning...'
   readcol,root+'_errs.txt',ell,errs,format='D'
   read_square_matrix,n_elements(ell),root+'_covar.txt',C 
   return
ENDIF

;MAD Read in simulations, get number of bins
readcol,infile,edge,cent,width,cl,run,format='D,D,D,D,I'
nbins=n_elements(where(run EQ 0))
ells=cent[where(run EQ 0)]

;MAD Build covarince matrix
C=dblarr(nbins,nbins)
FOR i=0L,max(run) DO BEGIN
  temp=cl[where(run EQ i)]
  C = C + (temp # temp)
ENDFOR
C=(1./(max(run)+1))*C

;MAD get diagonals for 1sigma errors
errs=dblarr(nbins)
FOR i=0L,n_elements(errs)-1 DO BEGIN
  FOR j=0L,n_elements(errs)-1 DO BEGIN
    IF (i EQ j) THEN errs[i]=sqrt(C[i,j])
  ENDFOR
ENDFOR

IF keyword_set(root) THEN BEGIN
   openw,1,root+'_covar.txt'
   formatstring=strarr(nbins)+'E,1x,'
   formatstring[0]='(E,1x,'
   formatstring[n_elements(formatstring)-1]='E)'
   formatstring=strjoin(formatstring)
   FOR i=0,nbins-1 DO BEGIN
      printf,1,C[i,*],format=formatstring
   ENDFOR
   close,1

   openw,1,root+'_errs.txt'
   FOR i=0L,n_elements(errs)-1 DO BEGIN
      printf,1,ells[i],errs[i]
   ENDFOR
   close,1
ENDIF

return
END
