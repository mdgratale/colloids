pro mgevplot, pos, eigenv, eigenn, scale=scale, vectorlength=vectorlength
;;; plot the eigen vectors at each particle in real space
;;; the ith row in the eigenv is the eigen vector of the ith eigen value
;scale=700.0
deigenv=eigenv(*,eigenn) ; the eigen value for the eigenn-th eigen value

npart=n_elements(pos[2,*])/2.0

x=pos[2,0:npart-1]
y=pos[2,npart:2*npart-1]

res=fltarr(4, npart)
	for i=0, npart-1 do begin
		res(0,i)=x(i)
		res(1,i)=y(i)  ;; x, y position
		res(2, i)=deigenv(i)
		res(3, i)=deigenv(i+npart)
	endfor

;	dev=stddev(res(2,*)^2+res(3,*)^2)
;	avg=mean(res(2,*)^2+res(3,*)^2)
;	w=where(res(2,*)^2+res(3,*)^2 ge avg )
;	temp=res(*, w)
temp=fltarr(5, npart)
temp(0:3,*)=res(0:3,*)
;temp(2,*)=(temp(2, *)-mean(temp(2,*))) ;; remove mode from drifting
;temp(3,*)=(temp(3, *)-mean(temp(2,*)))




veclength=sqrt(temp(2,*)^2+temp(3,*)^2)
avect=mean(veclength)
if keyword_set(vectorlength) then avect=vectorlength
;avect=median(veclength)
;temp(2,*)=(temp(2, *))/sqrt(veclength);; remove mode from drifting
;temp(3,*)=(temp(3, *))/sqrt(veclength) ;; scale the length

temp(2,*)=(temp(2, *))/(avect)
temp(3,*)=(temp(3, *))/(avect)
;plot_hist, veclength
;what
;maxv=max(veclength/avect)
;scale=50.0/maxv;; maximum as 50

temp(4, *)=atan(temp(3,*), temp(2,*))  ; angle
;stop
;window, 3,xsize=10*max(temp(0,*)), ysize=10*max(temp(1,*))
;plotp, temp
window,xsize=600,ysize=600
Z=max([x,y])
plot,x,y,xrange=[0,Z],yrange=[0,Z],psym=3
arrow, temp(0,*), temp(1,*), temp(0,*)+scale*temp(2,*), temp(1,*)+scale*temp(3,*), /data, hsize=-0.6, thick=2, /solid;, color=255
;legend, leg, /right, /top, charsize=1.5
;legend, leg2, /right, /center, charsize=
;return, temp
end


