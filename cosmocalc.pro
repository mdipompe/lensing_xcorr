;+
;  NAME:
;    cosmocalc
;  PURPOSE:
;    Get cosmological distances, etc.
;
;  USE:
;    res=cosmocalc(z,h=h,om=om,lambda=lambda)
;
;  INPUT:
;    z - redshift
;
;  OPTIONAL INPUT:
;    h - little h (H_0/100). Defaults to 0.71
;    om - omega_matter.  Defaults to 0.27
;    lambda - omega_lambda. Defaults to 0.73
;    
;  KEYWORDS:
;        
;  OUTPUT:
;    structure with 8 tags: d_h - the hubble distance                      
;                           d_m - the transverse comiving distance (Mpc)   
;                           d_L - the luminosity distance (Mpc)            
;                           d_a - the angular size distnace (Mpc)          
;                           d_c - the comoving distance (Mpc)              
;                           v_c - the comoving volume (Mpc^3)              
;                           t_l - the lookback time (Gyrs)                 
;                           t_h - the hubble time (Gyrs)                   
;  HISTORY:
;    2013 - Written - MAD (UWyo)
;    2015 - Cleaned and documented - MAD (UWyo)
;-
FUNCTION cosmocalc,z,h=h,om=om,lambda=lambda

;MAD Set speed of light, in km/s
c=2.99792458E5

;MAD Set defaults
IF ~keyword_set(h) THEN H0=71. ELSE H0=h*100.
IF ~keyword_set(om) THEN omega_m=0.27 ELSE omega_m=om
IF ~keyword_set(lambda) THEN omega_l=0.73 ELSE omega_l=lambda
omega_k=1.-omega_m-omega_l

;MAD Convert hubble constant to units of s^-1, calculate hubble distance
;(Mpc) and hubble time (yrs)
H0_conv=H0*(1./(3.0859E19))
d_h=(c/H0_conv)*(1./3.0859E19)
t_h=(1./H0_conv)*(1./3600.)*(1./24.)*(1./365.)

;MAD Integrate
IF (z LT 10.) THEN zvals=dindgen(z*100000.)/100000. ELSE zvals=dindgen(z*10000.)/10000.
E=1./SQRT(omega_m*((1+zvals)^(3.0))+Omega_k*((1+zvals)^(2.0))+omega_l)
E2=1./((1.+zvals)*SQRT(Omega_m*((1+zvals)^(3.0))+Omega_k*((1+zvals)^(2.0))+Omega_l))
y=int_tabulated(zvals,E,/double)
y2=int_tabulated(zvals,E2,/double)

;MAD Find comoving distance
d_c=d_h*y

;MAD Set d_m
IF Omega_k GT 0 THEN d_m=d_h*(1.0/SQRT(Omega_k))*SINH(SQRT(Omega_k)*(d_c/d_h))
IF Omega_k EQ 0 THEN d_m=d_c
IF Omega_k LT 0 THEN d_m=d_h*(1.0/SQRT((-1)*Omega_k))*SIN(SQRT((-1)*Omega_k)*(d_c/d_h))

;MAD Calculate angular size distance
d_a=d_m/(1.+z)

;MAD Calculate luminosity distance
d_L=d_m*(1.+z)

;MAD Calculate comoving volume within object
IF Omega_k GT 0 THEN v_c=((4.*!dpi*d_h^3.)/(2.*omega_k))*(((d_m/d_h)*SQRT(1.+omega_k*(d_m^2./d_h^2)))-((1./SQRT(abs(omega_k)))*ASINH(SQRT(abs(omega_k))*(d_m/d_h))))
IF Omega_k EQ 0 THEN v_c=((4.*!dpi)/3.)*d_m^3.
IF Omega_k LT 0 THEN v_c=((4.*!dpi*d_h^3.)/(2.*omega_k))*(((d_m/d_h)*SQRT(1.+omega_k*(d_m^2./d_h^2)))-((1./SQRT(abs(omega_k)))*ASIN(SQRT(abs(omega_k))*(d_m/d_h))))

;MAD Calculate lookback time
t_l=t_h*y2*(1e-9)

;MAD Build return structure
dists={d_h:0.D, d_m:0.D, d_L:0.D, d_a:0.D, d_c:0.D, v_c:0.D, t_l:0.D, t_h:0.D}
dists=replicate(dists,n_elements(z))
dists.d_h=d_h
dists.d_m=d_m
dists.d_L=d_L
dists.d_a=d_a
dists.d_c=d_c
dists.v_c=V_c
dists.t_l=t_l
dists.t_h=t_h*(1e-9)

return,dists
END


