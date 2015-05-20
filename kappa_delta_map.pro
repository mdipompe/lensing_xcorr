;+
;  NAME:
;    kappa_delta_map
;  PURPOSE:
;    Overlay portion of two HEALpix maps, one in color contours, one in line contours
;
;  USE:
;    kappa_delta_map,map1,map2,mask,lower_ra,upper_ra,lower_dec,upper_dec,'image.eps',split=split
;
;  INPUT:
;    map1 - HEALPix map to be plotted as color contours
;    map2 - HEALPix map to be plotted as line contours (same nside!)
;    mask - binary mask to apply to both maps
;    lra - lower limit of RA range (deg)
;    ura - upper limit of RA range (deg)
;    ldec - lower limit of DEC range (deg)
;    ura - upper limit of DEC range (deg)
;    imout - string of the eps image file name for output
;
;  OPTIONAL INPUTS:
;    split - value at which to switch linestyles for line contours of
;            map2.  Defaults to 0.
;
;  OUTPUT:
;    Writes eps file of image
;
;  HISTORY:
;    11-11-14 - Written - MAD (UWyo)
;-
PRO kappa_delta_map,map1,map2,mask,lra,ura,ldec,udec,imout,split=split

;MAD Define filled circle for plot
circsym

IF ~keyword_set(split) THEN split=0

;MAD deterine nside from number of pixels
nside = LONG(SQRT(n_elements(map1)/12))
pix=dindgen(12.*nside^2.)

;MAD Get pixel positions, convert to galactic coords
pix2ang_nest,nside,pix,theta,phi
undefine,pix
galra=phi*(180./!dpi)
galdec=90.-((180./!dpi)*theta)

;MAD convert to equatorial
euler,galra,galdec,ra,dec,2

;MAD free up some space...
undefine,galra
undefine,galdec
undefine,theta
undefine,phi

;MAD Set values to mean where data is masked
map1[where(mask EQ 0)]=mean(map1[where(mask EQ 1)])
map2[where(mask EQ 0)]=mean(map2[where(mask EQ 1)])

;MAD Limit arrays to your RA and DEC region
map1=map1[where((ra GE lra) AND (ra LE ura) AND (dec GE ldec) AND (dec LE udec))]
map2=map2[where((ra GE lra) AND (ra LE ura) AND (dec GE ldec) AND (dec LE udec))]
mask=mask[where((ra GE lra) AND (ra LE ura) AND (dec GE ldec) AND (dec LE udec))]
ra=ra[where((dec GE ldec) AND (dec LE udec))]
dec=dec[where((dec GE ldec) AND (dec LE udec))]
dec=dec[where((ra GE lra) AND (ra LE ura))]
ra=ra[where((ra GE lra) AND (ra LE ura))]


;Red/blue color bar
loadct,33

;MAD Use maximum number of levels for colors
n_levels=60
cols=fltarr(n_levels+1)
cols[0]=1.
FOR i=1,n_elements(cols)-1 DO BEGIN
 cols[i]=cols[i-1]+((300.-cols[0])/(n_levels+5))
ENDFOR
cols=round(cols)

;MAD Open file for writing
PS_start,filename=imout,/encapsul,xsize=17,ysize=13

;MAD Plot colored contours (kappa)
contour,map1,ra,dec,/irregular,$
    c_linestyle=0,/fill,nlevels=n_levels,c_colors=cols,$
    max_value=max(map1),min_value=min(map2),$
    xtit='RA',ytit='DEC',xra=[lra-1.,ura+1.],yra=[ldec-1,udec+6.],xsty=1,ysty=1,color=cgcolor('black'),$
    thick=2,xthick=6,ythick=6,charthick=3,xcharsize=2.,ycharsize=2.,position=[0.13,0.13,0.95,0.95]

;MAD Plot line contours (delta)
contour,map2[where(map2 LT split)],ra[where(map2 LT split)],dec[where(map2 LT split)],/irregular,$
        c_linestyle=2,/overplot,thick=3,$
        levels=[-0.5,-0.4,-0.3,-0.2,-0.1,0]
contour,map2[where(map2 GE split)],ra[where(map2 GE split)],dec[where(map2 GE split)],/irregular,$
        c_linestyle=0,/overplot,thick=3,$
        levels=[0,0.1,0.2,0.3,0.4,0.5]

;MAD Plot over bad data to make it white
oplot,ra[where(mask EQ 0)],dec[where(mask EQ 0)],psym=8,color=cgcolor('white'),symsize=0.2

;MAD Plot delta legend
loadct,0
delt=textoidl('\delta')
legend,[delt+' < '+strtrim(split,2),delt+' > '+strtrim(split,2)],linestyle=[2,0],$
       box=0,charsize=2,charthick=2,thick=3,position=[ura-35,udec+5.],/data

;MAD Plot kappa color bar
loadct,33
lowtickname=strtrim(round(min(map1)*100.)/100.,2)
hitickname=strtrim(round(max(map1)*100.)/100.,2)
lowtickname=string(lowtickname,format='(A5)')
hitickname=string(hitickname,format='(A4)')
kap=textoidl('\kappa')
colorbar,position=[0.2,0.91,0.5,0.93],$
         ticknames=[kap+'='+lowtickname,' ',' ',' ',' ',' ',kap+'='+hitickname],$
         annotatecolor='black',charsize=2.


PS_end

END
