PRO mgvecplot,xyz,vec,nt,mass1=mass1, scale=scale
;;; plot vectors, usually displacement, for number of time frames, also good for multiple species of particles
;;; xyz=positions, and masses if needed, can use track file
;;; vec=vector lengths, usually displacements
;;;mass1=mass cut off of species, nt= number of frames to image, scale=scale

if not keyword_set(nt) then nt=1.0
if not keyword_set(scale) then scale=1.0

ncol=n_elements(xyz[*,0])
xcol=0
ycol=1
bcol=2
np=n_elements(xyz[0,*])
maxtime=vec[0,*]
dx=vec[*,0:np-1]
dy=vec[*,np:2*np-1]
;stop
white=254.0+256.0*254.0+256.0*256.0*254.0
for ii=0,nt do begin
	if not keyword_set(mass1) then begin
		plot,xyz[xcol,*],xyz[ycol,*],psym=3
		arrow,xyz[xcol,*],xyz[ycol,*],xyz[xcol,*]+scale*dx[ii,*],xyz[ycol,*]+scale*dy[ii,*],color=white,/data,hsize=-0.65,thick=2.75
	endif else begin
		w1=where(xyz[bcol,*] gt mass1,n1,complement=w2,ncomplement=n2)
	    plot,xyz[xcol,*],xyz[ycol,*],psym=3
		arrow,xyz[xcol,w1],xyz[ycol,w1],xyz[xcol,w1]+scale*dx[ii,w1],xyz[ycol,w1]+scale*dy[ii,w1],color=white,/data,hsize=-0.7,thick=1.3,/solid
		arrow,xyz[xcol,w2],xyz[ycol,w2],xyz[xcol,w2]+0.4*scale*dx[ii,w2],xyz[ycol,w2]+0.4*scale*dy[ii,w2],color=254.0,/data,hsize=-0.7,thick=1.3,/solid
	endelse
	write_bmp,'displacements_t'+strcompress(long(ii),/remove_all)+'.bmp',tvrd(true=1)
endfor

END


