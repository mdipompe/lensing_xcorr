;+
;  NAME:
;    fit_crosspower
;  PURPOSE:
;    Fit a model to a CMB lensing cross-correlation
;
;  USE:
;    fit_crosspower,cl,ell,min_ell,max_ell,mod_ell,mod_cl,bias,err=err,covar=covar,plotout=plotout           
;
;  INPUT:
;    cl - cross-correlation coefficients
;    ell - ell values
;    min_ell - minimum ell value to fit
;    max_ell - maximum ell value to fit
;    mod_ell - model_ell values
;    mod_cl - model cross-power values
;
;  OPTIONAL INPUT:
;    err - 1-sigma error bars (set to very small if not provided)
;    covar  - Matrix of covariances (if a string, will read text file
;             containing covariance matrix).  If not set, just uses variance.
;    plotout - if provided, outputs hardcopy of plot. string name of file
;
;  KEYWORDS:
;    bz - if set assumes you input a model with b(z) and want to
;         rescale it (b_0).  Really just changes what is printed to screen
;
;  OUTPUT:
;    bias - the best fit bias
;    berr - error on the bias
;
;  OPTIONAL OUTPUT
;    minchi2 - the chi^2 of the best fit bias
;
;  HISTORY:
;    11-11-14 - Written - MAD (UWyo)
;    11-6-15 - Added minchi2 as output, b(z) option - MAD (Dartmouth)
;-

PRO fit_crosspower,cl,ell,min_l,max_l,mod_ell,mod_cl,bias,berr,$
                   err=err,covar=covar,plotout=plotout,$
                   minchi2=minchi2,bz=bz

IF (N_Params() LT 6) THEN message,'Syntax - fit_crosspower,cl,ell,min_l,max_l,''model.txt'',bias,[errors=errors,covar=covar,plotout=''plotfile.png'']'


;MAD get start time
st=systime(1)

;MAD Get filled circle for plot
circsym

;MAD Set errors to small value if not provided
IF ~keyword_set(err) then err=cl*0.000001

;MAD Read in covariance matrix if provided as file
IF keyword_set(covar) THEN BEGIN
   check=size(covar)
   IF ((check[0] EQ 0) AND (check[1] EQ 7)) THEN $
      read_square_matrix,n_elements(ell),covar,C ELSE C=covar
ENDIF

;MAD Limit to desired l range
xx=where((ell GE min_l) AND (ell LE max_l))
ell=ell[xx]
cl=cl[xx]
err=err[xx]
IF keyword_set(covar) THEN BEGIN
   C=C[*,xx]
   C=C[xx,*]
ENDIF

;MAD if no covariance provided, make a matrix with just
;variance elements
IF ~keyword_set(covar) THEN BEGIN
   C=dblarr(n_elements(ell),n_elements(ell))
   FOR i=0,n_elements(C[*,0])-1 DO BEGIN
      FOR j=0,n_elements(C[0,*])-1 DO BEGIN
         IF (i EQ j) THEN C[i,j]=err[i]^2.
      ENDFOR
   ENDFOR
ENDIF

;MAD Interpolate the model
newmod=interpol(mod_cl,mod_ell,ell)

;MAD Generate array of bias values for fit
b_guess=(findgen(10000.)*((5.-0.1)/10000.))+0.1

;MAD Initialize array of chisq vlaues for fit
chisq=dblarr(n_elements(b_guess))

;MAD Invert covariance matrix
C_inv=invert(C,/double)

;MAD Loop and calculate chisq
FOR i=0L,n_elements(b_guess)-1 DO BEGIN
   temp = ((cl-(newmod*b_guess[i]))) # C_inv # ((cl-(newmod*b_guess[i])))
   chisq[i]=chisq[i]+temp
ENDFOR

;MAD Find min chi^2, get best fit bias
minchi=where(chisq EQ min(chisq))
sol=b_guess[where(chisq EQ min(chisq))]

minchi2=chisq[minchi]

;MAD Find error on best fit bias, using delta chi^2 = 1
chi1=chisq[0:minchi]-min(chisq)
chi2=chisq[minchi:n_elements(chisq)-1]-min(chisq)
b_vals1=b_guess[0:minchi]
b_vals2=b_guess[minchi:n_elements(b_guess)-1]
xx=closest(chi1,1.)
b_err1=b_vals1[xx]
yy=closest(chi2,1.)
b_err2=b_vals2[yy]
b_err=(abs(sol[0]-b_err1)+abs(sol[0]-b_err2))*(0.5)
dchi_b=(chi1[xx]+chi2[yy])*0.5

;MAD Initialize plot hardcopy if needed
IF keyword_set(plotout) THEN PS_start,filename=plotout,xsize=11,ysize=9

;MAD Plot the results
ytit=textoidl('C_{l}^{\kappa q} \times 10^{6}')
xtit1=textoidl('l')
xtit2=textoidl('\theta (deg)')

plot,[0],[0],psym=0,$
     xtit=xtit1,ytit=ytit,xra=[6,2500],yra=[0.008,2],xsty=1,ysty=1,$
     thick=5,xthick=8,ythick=8,xticklen=0.000001,$
     charthick=1.0,xcharsize=1.8,ycharsize=1.8,/xlog,/ylog,$
     position=[0.19,0.17,0.93,0.87]

axis,6,0.008,xaxis=0,/xlog,xra=[6,2500],xsty=1,xticklen=0.02,$
     xcharsize=1.8,charthick=1.8,xthick=8

axis,(2.*!dpi/6.)*(180./!dpi),2.,xaxis=1,/xlog,xra=[(2.*!dpi/6.)*(180/!dpi),(2.*!dpi/2500.)*(180./!dpi)],$
  xsty=1,xtit=xtit2,xticklen=0.02,xcharsize=1.8,charthick=1.8,xthick=8

oplot,ell,cl*(1.e6),psym=8,color=cgcolor('dark grey'),symsize=1.5
oploterror,ell,cl*(1.e6),err*(1.e6),psym=3,color=cgcolor('dark grey')
oplot,mod_ell,mod_cl*sol[0]*1.e6,linestyle=1,thick=5


IF keyword_set(plotout) THEN PS_end,/png

IF ~keyword_set(bz) THEN BEGIN
   print,'Chi^2 of best fit bias is ' + strtrim(minchi2,2)
   print,'(error on b is based on delta chi of ',strtrim(dchi_b,2),')'
   print,'The best fit bias value is: '
   print,strtrim(sol,2)+' +/- '+strtrim(b_err,2)
   print,'(over the range l=',strtrim(min_l,2),'-',strtrim(max_l,2),')'
ENDIF ELSE BEGIN
   print,'Chi^2 of best fit b_0 is ' + strtrim(minchi2,2)
   print,'(error on b_0 is based on delta chi of ',strtrim(dchi_b,2),')'
   print,'The best fit b_0 value is: '
   print,strtrim(sol,2)+' +/- '+strtrim(b_err,2)
   print,'(over the range l=',strtrim(min_l,2),'-',strtrim(max_l,2),')'
ENDELSE

   
bias=sol
berr=b_err

return
END



