;+
;  NAME:
;    stack
;  PURPOSE:
;    Bin a HEALPix map, find locations within those bins of a
;    different map, and stack them.  Gives visual representation 
;    of how strongly maps are correlated.
;
;  USE:
;   stack,map1,map2,mask,stack1,stack2,outfile='plotfilename.eps',bins_edges=bin_edges,uperrs=uperrs,lowerrs=lowerrs
;
;  INPUT:
;    map1 - First HEALPix map to bin (usually a delta map)
;    map2 - HEALPix map to stack at locations within map1 bins
;    mask - HEALPix binary mask 
;
;  OPTIONAL INPUT
;    outplot - string name of output plot eps file
;    bin_edges - array of bin edges.  Defaults to bins in Geach et al. (2013)
;    uperrs - upper error bars in y direction    
;    lowerrs - lower error bars in y direction
;
;  OUTPUT:
;    stack1 - the mean value of map 1 in each bin
;    stack2 - the mean value of map 2 at the location of map 1 in each
;             bin
;    err1 - the standard deviation of map1 values in bins
;    plot
;
;  HISTORY:
;    11-16-14 - Written - MAD (UWyo)
;     2-20-17 - Minor improvements, added output for scatter in map1 bins
;-
PRO stack,map1,map2,mask,stack1,stack2,err1,outplot=outplot,bin_edges=bin_edges,uperrs=uperrs,lowerrs=lowerrs

;MAD Define filled circle for plot
circsym

;MAD Mask pixels
map1_good=map1[where(mask EQ 1)]
map2_good=map2[where(mask EQ 1)]

;MAD Bin map1, stack locations in map2.  Use standard deviations as errors
IF ~keyword_set(bin_edges) THEN $
   bin_edges=[-0.5,-0.27,-0.19,-0.12,-0.07,-0.01,0.04,0.1,0.17,0.27,0.5]

FOR i=0L,n_elements(bin_edges)-2 DO BEGIN
   xx=where((map1_good GT bin_edges[i]) AND (map1_good LT bin_edges[i+1]))
   IF (n_elements(stack2) EQ 0) THEN stack2=mean(map2_good[xx]) ELSE stack2=[stack2,mean(map2_good[xx])]
   IF (n_elements(signew1) EQ 0) THEN signew2=stddev(map2_good[xx]) ELSE signew2=[signew2,stddev(map2_good[xx])]

   IF (n_elements(stack1) EQ 0) THEN stack1=mean(map1_good[xx]) ELSE stack1=[stack1,mean(map1_good[xx])]
   IF (n_elements(signew1) EQ 0) THEN signew1=stddev(map1_good[xx]) ELSE signew1=[signew1,stddev(map1_good[xx])]
ENDFOR

err1=signew1

;MAD Make a plot
IF keyword_set(outplot) THEN BEGIN
   PS_start,filename=outplot,xsize=11,ysize=9,/encapsulate

   xtit=textoidl('\delta')
   ytit=textoidl('\kappa')

   plot,[0],[0],psym=0,/nodata,$
        xtit=xtit,ytit=ytit,xra=[-0.49,0.49],yra=[-0.014,0.02],xsty=1,ysty=1,$
        thick=5,xthick=8,ythick=8,$
        charthick=1.0,xcharsize=1.8,ycharsize=1.8,$
        position=[0.19,0.15,0.95,0.95],$
        ytickname=['-0.01',' ','0.00',' ','0.01',' ','0.02']
ENDIF ELSE BEGIN
   xtit=textoidl('\delta')
   ytit=textoidl('\kappa')

   plot,[0],[0],psym=0,/nodata,$
        xtit=xtit,ytit=ytit,xra=[-0.49,0.49],yra=[-0.014,0.02],xsty=1,ysty=1,$
        thick=3,xthick=4,ythick=4,$
        charthick=1.0,xcharsize=1.3,ycharsize=1.3,$
        position=[0.19,0.15,0.95,0.95],$
        ytickname=['-0.01',' ','0.00',' ','0.01',' ','0.02']
ENDELSE


oplot,stack1,stack2,psym=8,symsize=1.5,color=cgcolor('slate gray')
IF (keyword_set(uperr) AND keyword_set(lowerr)) THEN BEGIN
   oploterror,stack1,stack2,signew1,errneg*(-1.),/lobar,psym=3,color=cgcolor('slate gray')
   oploterror,stack1,stack2,signew1,errpos,/hibar,psym=3,color=cgcolor('slate gray')
ENDIF ELSE BEGIN
   oploterror,stack1,stack2,signew1,signew2,psym=3,color=cgcolor('slate gray')
ENDELSE

;MAD Overplot bin boundaries
y=[-10,10]
FOR i=0L,n_elements(bin_edges)-1 DO BEGIN
 IF keyword_set(outplot) THEN $
    oplot,[bin_edges[i],bin_edges[i]],y,linestyle=2,thick=8,color=cgcolor('gray') $
 ELSE oplot,[bin_edges[i],bin_edges[i]],y,linestyle=2,thick=4,color=cgcolor('gray')
ENDFOR

IF keyword_set(outplot) THEN PS_end

;MAD Test how well correlated the maps are
corr=r_correlate(stack1,stack2)
print,'STACK: Stacked maps correlation values r,p = ',corr

return
END
