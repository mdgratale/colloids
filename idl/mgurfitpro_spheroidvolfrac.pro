PRO mgurfitpro_spheroidvolfrac,ur,R2,L,phi,sigma,D,kappa,fitparams,fitvals,perr,chi,dof,pinfo=pinfo,scale=scale,error=error
;;fitparams=provided variable to hold fit parameters, fitvals=provided variable to hold fit function and data, perr=variable to hold fit parameter errors

if not keyword_set(scale) then scale=1.0

nr=n_elements(ur[0,*])

if not keyword_set(error) then error=fltarr(1,n_elements(ur[0,*]))+1.0

;rmax=max(ur[0,*])
;if rmax gt 2 then ur=ur[*,0:nr-2]
;if rmax gt 2 then error=error[0:n_elements(ur[0,*])-1]
;stop

X=ur[0,*]*scale
P=[R2,L,phi,sigma,D,kappa]

if not keyword_set(pinfo) then pinfo=replicate({fixed:0,limited:[0,0],limits:[0.D,0.D]},5)
fitparams=mpfitfun('mgurmodel_spheroidvolfrac',X,ur[1,*],error,P,parinfo=pinfo,perror=perr,yfit=fun,bestnorm=chi,dof=dof)

print,fitparams
print,perr
fitvals=[X,transpose(fun),ur[1,*],error]

END

