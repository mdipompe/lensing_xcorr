;+
;NAME
;   ransack
;
;PURPOSE
;   IDL wrapper for MANGLE function ransack, which populates
;   polygons with N random points
;
;USAGE
;   ransack,n,poly_file,outfile,seed=seed
;
;INPUTS
;   n - Number of random points to generate
;   polyfile - string name of file containing MANGLE polygons
;               to fill
;   outfile - string name of text file for ransack to output.
;             Contains two columns, RA and DEC
;
;OPTIONAL INPUTS
;   seed - Seed for random number generator.  Defaults to generating
;          random integer from 1 to 1000 using systime_seed.
;
;OUTPUTS
;
;HISTORY
;   1-5-15 - Written - MAD (UWyo)
;-
PRO ransack,n,polyfile,outfile,seed=seed

;MAD Generate seed if not set
IF ~keyword_set(seed) THEN BEGIN
   seed = ceil(randomu(systime_seed)*1000)
ENDIF

;MAD Make string command to run

cmd=[filepath('ransack', root_dir=getenv('IDLUTILS_DIR'), $
     subdir='bin'),'-c',strtrim(seed,2), '-r',strtrim(n,2),$
     polyfile,outfile]

spawn,cmd,/noshell

RETURN
END
