;+
;  NAME:
;    matter_power_spec
;
;  PURPOSE:
;    Generate an initialization file for CAMB that has many redshift
;    steps for the transfer function/matter power spectrum.  You
;    should first set up the initialization with all the necessary
;    CAMB parameters (cosmology, etc.).  This code will add lines for
;    the redshift steps and output files, and then run CAMB if keyword
;    is set.
;
;    Note that you should have CAMB installed, and the environment
;    variable CAMB_DIR set to the location of the camb executable.
;
;    Also note that CAMB will only run on up 150 redshifts at a time.
;    If you have more than this you should split up your runs.  This
;    procedure will quit if you have more than 150.
;
;  USE:
;    matter_power_spec,'paramfile','out_paramfile',z_step,min_z,max_z,/camb
;
;  INPUT:
;    paramfile - the CAMB initialization file with general parameters
;                set
;    outparamfile - the name for the output parameter file with all of
;                   the necessary zs added in
;    zstep - the step size in z
;    zmin - the smallest value of z
;    zmax - the largest value of z
;
;  KEYWORDS:
;    camb - if set, will run CAMB on the generated parameter file
;
;  OUTPUT:
;    Writes new parameter file with general settings preserved, but
;    lines for each transfer function redshift.  If /camb is set, will
;    produce CAMB output files.
;
;  HISTORY:
;    5-20-15 - Written - MAD (UWyo)
;-
PRO matter_power_spec,paramfile,outparamfile,zstep,minz,maxz,camb=camb

;MAD Generate z array (reverse is because CAMB requires calculations
;MAD to be done in decreasing z)
z=(findgen((maxz-minz)/zstep)*zstep)+minz
z=[z,max(z)+zstep]
z=reverse(z)


;MAD Check number of redshifts
IF (n_elements(z) GT 150) THEN BEGIN
   print,'**** CAMB will only run on up to 150 redshifts, you have ' + strtrim(n_elements(z),2) + ' ****'
   print,'Exiting...'
   stop
ENDIF

;MAD generate string array of lines from initial param file, find
;MAD locations of places to add lines later
openr,1,paramfile
line=''
i=0
WHILE ~EOF(1) DO BEGIN
   readf,1,line
   subline=strsplit(line,'=',/extract)
   IF (strtrim(subline[0],2) EQ 'transfer_num_redshifts') THEN nline=i
   IF (strtrim(subline[0],2) EQ 'transfer_redshift(1)') THEN zline=i
   IF (strtrim(subline[0],2) EQ 'transfer_filename(1)') THEN tfileline=i
   IF (strtrim(subline[0],2) EQ 'transfer_matterpower(1)') THEN mpsline=i
   IF (n_elements(params) EQ 0) THEN params=line ELSE $
      params=[params,line]
   i=i+1
ENDWHILE
close,1


;MAD Start generating new param file
;MAD Write general lines up to number of zs for transfer function
openw,1,outparamfile
FOR i=0L,nline-1 DO BEGIN
   printf,1,params[i],format='(A)'
ENDFOR

;MAD Write number of zs
printf,1,'transfer_num_redshifts = ' + strtrim(n_elements(z),2),format='(A)'
printf,1,params[nline+1],format='(A)'

;MAD Write transfer function zs
FOR i=0L,n_elements(z)-1 DO BEGIN
   string='transfer_redshift(' + strtrim(i+1,2) + $
          ')    = ' + strtrim(z[i],2)
   printf,1,string,format='(A)'
ENDFOR

;MAD Write transfer function output file names
FOR i=0L,n_elements(z)-1 DO BEGIN
   string='transfer_filename(' + strtrim(i+1,2) + $
          ')    = ' + 'transfer_out' + $
          strtrim(z[i],2) + '.dat'
   printf,1,string,format='(A)'
ENDFOR

printf,1,params[tfileline+1],format='(A)'

;MAD Write Matter Power Spectrum output file names
FOR i=0L,n_elements(z)-1 DO BEGIN
   string='transfer_matterpower(' + strtrim(i+1,2) + $
          ')    = ' + 'matterpower' + $
          strtrim(z[i],2) + '.dat'
   printf,1,string,format='(A)'
ENDFOR

;MAD Write rest of general lines
FOR i=mpsline+1,n_elements(params)-1 DO BEGIN
   printf,1,params[i],format='(A)'
ENDFOR
close,1


;MAD Run CAMB, if desired
IF keyword_set(camb) THEN BEGIN
   ;MAD Copy over file needed to run CAMB in current directory
   cmd=['cp',filepath('HighLExtrapTemplate_lenspotentialCls.dat', root_dir=getenv('CAMB_DIR')),'.']
   spawn,cmd,/noshell

   ;MAD Run CAMB with new ini file
   cmd=[filepath('camb', root_dir=getenv('CAMB_DIR')),strtrim(outparamfile,2)]
   print,cmd
   spawn,cmd,/noshell
ENDIF


END
