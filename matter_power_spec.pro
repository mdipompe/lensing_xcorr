;+
;  NAME:
;    matter_power_spec
;
;  PURPOSE:
;    Generate an initialization file for CAMB that has many redshift
;    steps for the transfer function/matter power spectrum.  Then will
;    use CAMB4IDL to get the matter power spectrum at each z.  You
;    should set whatever defaults you want in your paramfile, but
;    cosmology and redshifts will be reset as supplied here.
;
;    Note that you should have CAMB installed, and the environment
;    variable CAMB_DIR set to the location of the camb executable.
;
;    Also note that CAMB will only run on up 150 redshifts at a time.
;    If you have more than this you should split up your runs.  This
;    procedure will quit if you have more than 150.
;
;  USE:
;    matter_power_spec,'paramfile',zarray,cosmology=....
;
;  INPUT:
;    paramfile - the CAMB initialization file with general parameters
;                set
;    zarray - array of z values to get power spectrum at
;
;  KEYWORDS:
;    h0 - little h (default 0.702)
;    omega_b - baryon fraction (default 0.046)
;    omega_dm - CDM fraction (default 0.229)
;    omega_l - lambda (default 0.725)
;    maxk - max wavenumber to run to (defaults to CAMB default of 2)
;
;  OUTPUT:
;    CAMB parameter files and CAMB output
;
;  HISTORY:
;    5-20-15 - Written - MAD (UWyo)
;    7-16-15 - Modified to use CAMB4IDL - MAD (UWyo)
;-
PRO matter_power_spec,paramfile,zarray,h0=h0,omega_b=omega_b,omega_dm=omega_dm,omega_l=omega_l,maxk=maxk

;MAD Set default cosmology
IF ~keyword_set(h0) THEN h0=0.702
IF ~keyword_set(omega_b) THEN omega_b=0.046
IF ~keyword_set(omega_dm) THEN omega_dm=0.225
IF ~keyword_set(omega_l) THEN omega_l=0.725

IF ~keyword_set(maxk) THEN maxk=2

;MAD Reverse zs if needed
IF (zarray[n_elements(zarray)-1] GT zarray[0]) THEN $
   revz=reverse(zarray) ELSE revz=zarray


;MAD Check number of redshifts
IF (n_elements(revz) GT 150) THEN BEGIN
   print,'**** CAMB will only run on up to 150 redshifts, you have ' + strtrim(n_elements(revz),2) + ' ****'
   print,'Exiting...'
   stop
ENDIF

tfile=strarr(n_elements(revz))
mpfile=strarr(n_elements(revz))
FOR j=0L,n_elements(revz)-1 DO BEGIN
   tfile[j]='transfer_' + strtrim(revz[j],2) + '.dat'
   mpfile[j]='matterpower_' + strtrim(revz[j],2) + '.dat'
ENDFOR

res=camb4idl(/runcamb, paramfile=paramfile,output_root='camb', $
             get_scalar='T',get_transfer='T', $
             get_tensor='F',do_lensing='T', $
             transfer_kmax=maxk, $
             transfer_num_redshifts=n_elements(revz), $
             transfer_redshift=revz, $
             transfer_filename=tfile, transfer_matterpower=mpfile, $
             use_physical='F', $
             hubble=h0*100., $
             omega_baryon=omega_b, $
             omega_cdm=omega_dm, $
             omega_lambda=omega_l)



END
