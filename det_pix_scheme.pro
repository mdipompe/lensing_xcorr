;+
;  NAME:
;    det_pix_scheme
;
;  PURPOSE:
;    Determine the pixelization scheme of a Mangle polygon file
;
;  USE:
;    det_pix_scheme,'polygonfile.ply'
;
;  INPUT:
;    polyfile - string name of the mangle polygon file
;             
;  OUTPUT:
;    Returns string with the pixelization scheme (or 'none' if not pixelized)
;
;  HISTORY:
;    5-20-15 - Written - MAD (UWyo)
;-
FUNCTION det_pix_scheme,polyfile

;MAD Open polygon file, loop through until pixelization scheme is
;found or get to the polygons (and no pixelization is present)
openr,1,polyfile
txt=''
flag=0
WHILE (flag EQ 0) DO BEGIN
   readf,1,txt
   xx=strsplit(txt,' ',/extract)
   IF (xx[0] EQ 'pixelization') THEN BEGIN
      scheme=xx[1]
      flag=1
   ENDIF
   IF (xx[0] EQ 'polygon') THEN flag=1
ENDWHILE

;MAD If pixelization scheme not found, return to 'none'
IF (n_elements(scheme) EQ 0) THEN scheme='none'
   
return,scheme
END
