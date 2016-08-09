function mgdisplaceme,otra,lifetime=lifetime,dim=dim,res=res

if not keyword_set(dim) then dim=2
;if not keyword_set(res) then res=1.0

tra=otra
ncol=n_elements(tra(*,0))
xcol=0
ycol=1
bcol=2   ;; brightness column if needed
idcol=ncol-1
timecol=ncol-2
npart=max(tra(idcol, *))+1 ; number of particles
maxtime=max(tra(timecol,*))+1 ; longest time
if not keyword_set(lifetime) then lifetime=round(0.9*maxtime)

traend=n_elements(tra(idcol,*)) ; the last index of the trackfile plus one
ididx=tra(idcol,*)
ididx=ididx-shift(ididx,1)
wp=where(ididx(*) ne 0.0,nwp)
wp=[wp, traend]  ;; the ith particle is between wp(i) and wp(i+1)-1

counter=0
temp=fltarr(nwp*dim, maxtime) ; twice the x size for both x and y directions--memory intensive!!
table=fltarr(dim+3,nwp)      ;; position table if needed, [id,x,y,b,life], changed by MG
weight=1.0
xcount=0
for i=0, nwp-1 do begin
	life=wp(i+1)-wp(i)-1
	if life ge lifetime then begin
	      ;xpos=fltarr(1,maxtime)
        time=tra(timecol,wp(i):wp(i+1)-1)
        ;xpos(round(time))=tra(xcol,wp(i):wp(i+1)-1)
        xavg=mean(tra(xcol,wp(i):wp(i+1)-1))
    	  if not keyword_set(res) then begin
    	     temp[counter,round(time)]=tra(xcol,wp(i):wp(i+1)-1)-xavg
    	  endif else begin
    	     temp[counter,round(time)]=res*round((tra(xcol,wp(i):wp(i+1)-1)-xavg)/res)
        endelse
        table(0,xcount)=tra(idcol,wp(i)) ;track ID
    	  table(1,xcount)=xavg
    	  table(3,xcount)=mean(tra(bcol,wp(i):wp(i+1)-1))
    	  table(4,xcount)=life
	      xcount=xcount+1
      	counter=counter+1
	endif
endfor
;stop
ycount=0
for i=0, nwp-1 do begin
	life=wp(i+1)-wp(i)-1
	if life ge lifetime then begin
  	  ;ypos=fltarr(1,maxtime)
      time=tra(timecol,wp(i):wp(i+1)-1)
	    ;ypos(round(time))=tra(ycol,wp(i):wp(i+1)-1)
      yavg=mean(tra(ycol,wp(i):wp(i+1)-1))
      if not keyword_set(res) then begin
           temp[counter,round(time)]=tra(ycol,wp(i):wp(i+1)-1)-yavg
      endif else begin
           temp[counter,round(time)]=res*round((tra(ycol,wp(i):wp(i+1)-1)-yavg)/res)
      endelse
    	;table(0,ycount)=tra(idcol,wp(i)) ;track ID
  	  table(2,ycount)=yavg
    	counter=counter+1
	    ycount=ycount+1
	endif
endfor

if dim eq 3 then begin
for i=0, nwp-1 do begin
  	time=tra(timecol, wp(i):wp(i+1)-1)
    	bavg=mean(long(tra(bcol,wp(i):wp(i+1)-1)))
  	temp[counter,round(time)]=tra(bcol,wp(i):wp(i+1)-1)-bavg
    	table(0,counter)=tra(idcol,wp(i)) ;track ID
  	table(1,counter)=bavg
    	counter=counter+1
endfor
endif

if xcount ne ycount then print,'Houston, we have a problem.'

;table2=fltarr(5,xcount)
;table2[0:2,*]=table[0:2,0:xcount-1]

dmx=temp[0:counter-1,*] ;; displacment matrix
write_gdf, dmx, 'displacements_big.gdf'

;dx=dmx[0:xcount-1,*]
;dy=dmx[xcount:2*xcount-1,*]
;dis=sqrt(dx^2+dy^2)
;for jj=0,xcount-1 do begin
;	m=max(dis[jj],mid)
;	table2[3,jj]=dx[mid]
;	table2[4,jj]=dy[mid]
;endfor

write_gdf,table[*,0:xcount-1],'avg_pos_big.gdf'

return,dmx
end
