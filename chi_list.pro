;+
;  NAME:
;    chi_list
;  PURPOSE:
;    Make an array of redshifts and their corresponding co-moving distances
;
;  USE:
;    chi_list,maxz,z,chi,zstep=zstep,omega_m=omega_m,omega_l=omega_l,h0=h0,outfile='z_chi_list.txt'
;
;  INPUT:
;    maxz - the largest redshift value
;    
;  OPTIONAL INPUTS:
;    zstep - the z stepsize (defaults to 0.01)
;    omega_m - Omega_matter, defaults to 0.273
;    omega_l - Omega_lambda, defualts to 0.727
;    h - little h (H0/100), defaults to 0.702
;
;  OUTPUT:
;    z - array of redshifts
;    chi - array of comoving distances
;    outfile - if supplied, writes text file of z, chi values
;
;  HISTORY:
;    11-11-14 - Written - MAD (UWyo)
;-
PRO chi_list,maxz,z,chi,zstep=zstep,omega_m=omega_m,omega_l=omega_l,h0=h0,outfile=outfile

;MAD If output file already exists, don't run just read it in
IF ~keyword_set(outfile) THEN check='' ELSE check=file_search(outfile)
IF (check NE '') THEN BEGIN
   print,'CHI_LIST: Output file already exists, reading in and returning...'
   readcol,outfile,z,chi
   return
ENDIF

;MAD Set defaults
IF ~keyword_set(omega_m) THEN omega_m=0.275
IF ~keyword_set(omega_l) THEN omega_l=0.725
IF ~keyword_set(h0) THEN h0=0.702
IF ~keyword_set(zstep) THEN zstep=0.01

;MAD Generate array of z values
z=(findgen(maxz/zstep)/(1./zstep))+zstep
IF (max(z) LT maxz) THEN z=[z,max(z)+zstep]

;MAD Call cosmology calculator to get chi for each z
;chi=cosmo_gen('c',omega_m,omega_l,h0,z)  ;Old and slower (but slightly more accurate...)
chi=fltarr(n_elements(z))
FOR i=0L,n_elements(z)-1 DO BEGIN
   d=cosmocalc(z[i],h=h0,om=omega_m,lambda=omega_l)
   chi[i]=d.d_c
ENDFOR

;MAD Write file if needed
IF keyword_set(outfile) THEN BEGIN
   openw,1,outfile
   printf,1,';z     chi   (omega_m='+strtrim(omega_m,2)+$
          '  omega_l='+strtrim(omega_l,2)+'  h0='+strtrim(h0,2)
   FOR i=0L,n_elements(chi)-1 DO BEGIN
      printf,1,z[i],chi[i]
   ENDFOR
   close,1
ENDIF

return
END
