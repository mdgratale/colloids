;02-11-11, Matt Gratale
;simulation program for crystal cluster DOS

FUNCTION simcrystalclusterDOS,L,D,k,quiet=quiet
;L=length of lattice,D=diameter of lattice, k=spring value

H=max([L,D])
H=float(H)

np=floor(1.134*H^2)  ;number of particles
t0=fltarr(4,np) ;array to hold x,y,k,id
;set x-postions
for ijk=0,np-1 do begin
  t0[0,ijk]=ijk mod H
endfor
;set y-positions
for j=0,(np/H)-1 do begin
  t0[1,j*H:(j+1)*H-1]=j*0.866
  if j mod 2 eq 1 then t0[0,j*H:(j+1)*H-1]=t0[0,j*H:(j+1)*H-1]+0.5
endfor

t=eclip(t0,[0,0,D],[1,0,L]) ;clip to desired width and length
N=n_elements(t[0,*])

t[2,*]=k

t[3,*]=findgen(1,N) ;give ID numbers

table=fltarr(3,2*N)
table[0,0:N-1]=t[3,*] ;ID
table[0,N:2*N-1]=t[3,*] ;ID
table[1,0:N-1]=t[2,*] ;k
table[1,N:2*N-1]=t[2,*] ;k
table[2,0:N-1]=t[0,*] ;x-pos.
table[2,N:2*N-1]=t[1,*] ;y-pos.

write_gdf,table,'positions'+strcompress(L,/remove_all)+'x'+strcompress(D,/remove_all)+'k'+strcompress(k,/remove_all)+'.gdf'

x=fltarr(N,N)
y=fltarr(N,N)
for ix=0,N-1 do x[ix,*]=t[0,ix] ;x-coordinates
for jy=0,N-1 do y[jy,*]=t[1,jy] ;y-coordinates
dx=x-transpose(x) ;delta-x
dy=y-transpose(y) ;delta-y
dist=sqrt(dx^2+dy^2)           ;; distance between particles
;kmatrix=(randomu(seed,2*N,2*N)-0.5)*2*0.0001 ;kmatrix array
kmatrix=fltarr(2*N,2*N)
;find neighbors and set springs
for ii=0,N-1 do begin
	nbr=where(dist[ii,*] gt 0 and dist[ii,*] le 1.0,nn)
	for jj=0,nn-1 do begin
	rxd=dx[ii,nbr[jj]]/dist[ii,nbr[jj]]
	ryd=dy[ii,nbr[jj]]/dist[ii,nbr[jj]]
	;first kxx
	;kmatrix[ii,nbr[jj]]=-0.5*(t[2,ii]+t[2,nbr[jj]])*(1+(randomn(seed)/sqrt(N)))*rxd^2
	kmatrix[ii,nbr[jj]]=-k*rxd^2
	;now kyy
	;kmatrix[ii+N,nbr[jj]+N]=-0.5*(t[2,ii]+t[2,nbr[jj]])*(1+(randomn(seed)/sqrt(N)))*ryd^2
	kmatrix[ii+N,nbr[jj]+N]=-k*ryd^2
	;now kxy
	;kmatrix[ii+N,nbr[jj]]=-0.5*(t[2,ii]+t[2,nbr[jj]])*(1+(randomn(seed)/sqrt(N)))*rxd*ryd
	kmatrix[ii+N,nbr[jj]]=-k*rxd*ryd
	;now kyx
	;kmatrix[ii,nbr[jj]+N]=-0.5*(t[2,ii]+t[2,nbr[jj]])*(1+(randomn(seed)/sqrt(N)))*rxd*ryd
	kmatrix[ii,nbr[jj]+N]=-k*rxd*ryd
	endfor
;stop
endfor
;stop
kmatrix=(kmatrix+transpose(kmatrix))/2.0 ;make sure symmetric
;make sure rows add to zero
for bb=0,2*N-1 do begin
	tot=total(kmatrix[bb,*])
	kmatrix[bb,bb]=kmatrix[bb,bb]-tot
endfor

kmatrix=(kmatrix+transpose(kmatrix))/2.0 ;make sure symmetric
write_gdf,kmatrix,'kmatrix'+strcompress(L,/remove_all)+'x'+strcompress(D,/remove_all)+'k'+strcompress(k,/remove_all)+'.gdf'

;eigen=la_eigenql(kmatrix,eigenvectors=egv,failed=fucked)
eigen=eigenql(kmatrix,eigenvectors=egv)
write_gdf,eigen,'eigenvalues'+strcompress(L,/remove_all)+'x'+strcompress(D,/remove_all)+'k'+strcompress(k,/remove_all)+'.gdf'
write_gdf,egv,'eigenvectors'+strcompress(L,/remove_all)+'x'+strcompress(D,/remove_all)+'k'+strcompress(k,/remove_all)+'.gdf'
;write_gdf,fucked,'failed_eigenvalues.gdf'

return,eigen
End

