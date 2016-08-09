;05-23-10, Matt Gratale
;function to display xyz modes

pro display_modes_xymg,xy,ev,rad,i,nt
;xy=average positions given by DOS code, ev=eigenvectors, rad=desired radius of circles, i=mode
;nt=number of frames made, if movie desired

if not keyword_set(nt) then nt=1

	for t=0,nt-1 do begin

np=n_elements(ev(*,0))/2.0

;pf=ev
;for j=0, 2*np-1 do begin
;	pf(*,j)=ev(*,j)/sqrt(total(ev(*,j)*ev(*,j)))
;endfor

evp=fltarr(6,np)
evp(0,*)=xy(2,0:np-1) ;avg. x
evp(1,*)=xy(2,np:2*np-1) ;avg. y
evp(2,*)=xy(1,0:np-1) ;mass
evp(3,*)=0 ;column of zeros for plotellipse (plotting cirlces here)
evp(5,*)=rad

white=254.0+256.0*256.0*254.0+256.0*254.0
ps=where(evp[2,*] gt 1.0,complement=nipa,numps)
colo=fltarr(1,np)
if numps gt 0 then begin
	colo(ps)=254.0+256.0*0.0+256.0*256.0*0.0
	colo(nipa)=254.0+256.0*254.0+256.0*256.0*150.0
endif else begin
	colo(*)=254.0+256.0*254.0+256.0*256.0*150.0
endelse

plotellipse,evp(0:2,*),a=evp(5,*),b=evp(5,*),colo=colo,scal=rad
evplot, evp(0,*),evp(1,*), ev, i, scale=25*cos(2*!pi*t/nt),vectorlength=0.1
invfreq=long(2*np-1)
write_bmp,'20000frames_mode'+strcompress(invfreq-long(i),/remove_all)+'_'+strcompress(long(t),/remove_all)+'.bmp',tvrd(true=1)

	endfor

end
