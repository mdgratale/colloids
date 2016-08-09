;02-11-11, Matt Gratale
;simulation program for crystal DOS

function simcrystalDOS,L,k1,k2,frac,numit,quiet=quiet
;L=length of lattice, k1,k2=spring values, frac=fraction of particles with k2, numit=number of iterations

if not keyword_set(numit) then numit=1

N=floor(1.134*L^2)  ;number of particles
t=fltarr(4,N) ;array to hold x,y,'mass',id
;set x-postions
for ijk=0,N-1 do begin
  t[0,ijk]=ijk mod L
endfor
;set y-positions
for j=0,(N/L)-1 do begin
  t[1,j*L:(j+1)*L-1]=j*0.866
  if j mod 2 eq 1 then t[0,j*L:(j+1)*L-1]=t[0,j*L:(j+1)*L-1]+0.5  
endfor
;egtemp=fltarr(numit,2*N) ;array to hold eigenvalues
pftemp=fltarr(3,2*N) ;array to hold participation fractions

for i=0,numit-1 do begin

seed=i

t[2,*]=randomu(seed,N) ;randomize mass
if frac gt 0.0 and frac lt 1.0 then begin
w=where(t[2,*] le frac,complement=z,nw)
t[2,w]=k2
t[2,z]=k1
endif
if frac eq 0.0 then begin
t[2,*]=k1
endif
if frac eq 1.0 then begin
t[2,*]=k2
endif
t[3,*]=findgen(1,N) ;give ID numbers

table=fltarr(3,2*N)
table[0,0:N-1]=t[3,*] ;ID
table[0,N:2*N-1]=t[3,*] ;ID
table[1,0:N-1]=t[2,*] ;'mass'
table[1,N:2*N-1]=t[2,*] ;'mass'
table[2,0:N-1]=t[0,*] ;x-pos.
table[2,N:2*N-1]=t[1,*] ;y-pos.
w2=where(table[1,*] eq k1,complement=z2,nw2)
if nw2 ne 0 then begin
	table[1,w2]=1.0
endif
if nw2 lt 2*N then begin
	table[1,z2]=1.1
endif
;write_gdf,table,'positions.gdf'

x=fltarr(N,N)
y=fltarr(N,N)
for ix=0,N-1 do x[ix,*]=t[0,ix] ;x-coordinates
for jy=0,N-1 do y[jy,*]=t[1,jy] ;y-coordinates
dx=x-transpose(x) ;delta-x
dy=y-transpose(y) ;delta-y
dist=sqrt(dx^2+dy^2)           ;; distance between particles
kmatrix=(randomu(seed,2*N,2*N)-0.5)*2*0.0001 ;kmatrix array
;find neighbors and set springs
for ii=0,N-1 do begin
	nbr=where(dist[ii,*] gt 0 and dist[ii,*] le 1.0,nn)
	for jj=0,nn-1 do begin
	rxd=dx[ii,nbr[jj]]/dist[ii,nbr[jj]]
	ryd=dy[ii,nbr[jj]]/dist[ii,nbr[jj]]
	;first kxx
	kmatrix[ii,nbr[jj]]=-0.5*(t[2,ii]+t[2,nbr[jj]])*(1+(randomn(seed)/sqrt(N)))*rxd^2
	;now kyy
	kmatrix[ii+N,nbr[jj]+N]=-0.5*(t[2,ii]+t[2,nbr[jj]])*(1+(randomn(seed)/sqrt(N)))*ryd^2
	;now kxy
	kmatrix[ii+N,nbr[jj]]=-0.5*(t[2,ii]+t[2,nbr[jj]])*(1+(randomn(seed)/sqrt(N)))*rxd*ryd
	;now kyx
	kmatrix[ii,nbr[jj]+N]=-0.5*(t[2,ii]+t[2,nbr[jj]])*(1+(randomn(seed)/sqrt(N)))*rxd*ryd
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
;write_gdf,kmatrix,'kmatrix.gdf'

;eigen=la_eigenql(kmatrix,eigenvectors=egv,failed=fucked)
eigen=eigenql(kmatrix,eigenvectors=egv)
;write_gdf,eigen,'eigenvalues.gdf'
;write_gdf,egv,'eigenvectors.gdf'
;write_gdf,fucked,'failed_eigenvalues.gdf'

;egtemp[i,*]=eigen
pf=sim_partfrac(eigen,egv,table)
pftemp=pftemp+pf
if not keyword_set(quiet) then print,i+1
endfor

;egavg=egtemp/numit ;find averages
pfavg=pftemp/numit

;write_gdf,egavg,'avg_eigenvalues.gdf'
;write_gdf,pfavg,'avg_partfrac.gdf'
;write_text,pfavg,'50-50_k1-50_k2-90_partfrac.txt'

return,pfavg
end

