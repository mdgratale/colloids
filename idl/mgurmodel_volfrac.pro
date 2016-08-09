;;function g(r) for mpfitfun

FUNCTION mgUrmodel_volfrac,r,param; param=[B,L,phi,sigma,D]
;;want things in units of nanometers
;;r=distances, L=micelle length in nanometers, nb=micelle number fraction, sigma=polydispersity

B=param[0]
;B=1640.0
L=param[1]
phi=param[2]
sigma=param[3]
D=param[4] ;width of micelles
D2=D*D
;B=1570. ;colloid radius in nanameters
;phi=nb*!pi*L/4
;sigma=B*poly
;sigma=50.0
nr=n_elements(r)
dr=r[1]-r[0]
;t=((findgen(1,75)/74.0)*8*sigma)-(4*sigma) ;set range of values for polydispersity
nt=1001.0
dG=dr/10.0
t=(findgen(1,nt)*dG)
t=t-t[(nt-1)/2.0]
p=exp(-(t)^2/2/sigma/sigma)/sqrt(2*!pi)/sigma ;gaussian kernel
p=p/total(p) ;rescale kernel to sum of 1
dia=t+B ;find particle diameters
urpoly=fltarr(nt,nr) ;array to hold particle polydispersities
grpoly=urpoly

for ii=0,nt-1 do begin
	w0=where(r ge dia(ii) and r le dia(ii)+L,n0) ;find h values within potential well
	winf=where(r lt dia(ii),ninf)
	ugr=fltarr(3,nr) ;array to hold resuts
	ugr[0,*]=r
	U=0
	if n0 gt 0 then begin
		h=r[w0]-dia(ii)
		U=-1.0*(2./3.)*phi*(0.5*dia(ii)*L/D2)*((1-(h/L))^3) ;calculate U(r)
		ugr[2,w0]=U
	endif
	if ninf gt 0 then ugr[2,winf]=1.0/0.0 ;set "overlap" r values to infinity
	ugr[1,*]=exp(-1.0*ugr[2,*]) ;invert to g(r)
	;ugr[1,winf]=0
	;if ii eq 37 then stop
	urpoly[ii,*]=ugr[2,*]*p[ii]
	grpoly[ii,*]=ugr[1,*]*p[ii] ;list of g(r)*probability
endfor

smoogr=total(grpoly,1)
smooUr=total(urpoly,1)

;smooUr=-1.0*alog(smoogr)
;fin=finite(smooUr)
;z=where(fin eq 0)
;smooUr[z]=1e4
;smoogr2=convol(ugr[1,*],p,/edge_truncate,/normalize)
;smoogr=ugr[0:1,*]
;smoogr[1,*]=gaus_convol(ugr[0,*],ugr[1,*],sigma,nsigma=3)

;smooUr=ugr[0:1,*]
;smooUr[1,*]=gaus_convol(ugr[0,*],ugr[2,*],sigma,nsigma=10)
;stop
;ws=where(smoogr[0,*] lt 1700); and smoogr[1,*] lt 0.5,ns)
;smoogr[1,ws]=1.0
;stop

return,-1.0*alog(smoogr)
END