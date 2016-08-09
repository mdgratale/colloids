;;; density of states
;;; packaged programs using matrices

function DOSphmg, otra, lifetime, bri=bri, mass1=mass1;, res=res;,ev=ev
;; convert track to a displacement matrix, not weighted
;; calculate the correlation between displacements
;; calculate the eigen values from correlation matrix
;; otra is the track file, lifetime is the number of frames the particle needs to be
;; to be included in the calculation (short lifetime messes up correlation, choose at least 90% of your number of frames),
;; sm is the smooth factor in the remove motion (drfiting creats false correlations)
;; insufficient statistics may result in the failure of solving the matrix
;; if only the modes are needed run this function directly
;; if you need some intermediate results, positions, matrix, vector... uncomment certain lines and save them
;;  I kept the weighted lines commented, as it is a issue remains to be worked out theoretically, though the results should not be quite different


;;!!!Many changes made by Matt Gratale, originally written to find density of states of NIPA-PS system, don't say you weren't warned
;;bri is brightness cutoff to seperate different species/masses of particles
;;mass1 is the mass of particles with brightness exceeding bri
;;res is assumed resolution

if not keyword_set(bri) then bri=0.0
if not keyword_set(mass1) then mass1=1.0
;if not keyword_set(res) then res=1.0

;mot=motion(otra)
;tra=rm_motion(otra, mot, smooth=sm) ;; remove drifting
;otra=0                              ;; clear original track
 tra=otra
 ncol=n_elements(tra(*,0))
 xcol=0
 ycol=1
 bcol=2   ;; brightness column if needed
 idcol=ncol-1
 timecol=ncol-2
npart=max(tra(idcol, *))+1 ; number of particles
maxtime=max(tra(timecol,*))+1 ; longest time

traend=n_elements(tra(idcol,*)) ; the last index of the trackfile plus one
ididx=tra(idcol,*)
ididx=ididx-shift(ididx,1)
wp=where(ididx(*) ne 0.0,nwp)
wp=[wp, traend]  ;; the ith particle is between wp(i) and wp(i+1)-1

counter=0
temp=fltarr(npart*2, maxtime) ; twice the x size for both x and y directions--memory intensive!!
table=fltarr(3, npart*2)      ;; position table if needed, x0...xn, y0...yn
weight=1.0

for i=0, nwp-1 do begin;npart-1 do begin change by PJY
	life=wp(i+1)-wp(i)-1
	if life ge lifetime then begin
        if tra(bcol, wp(i)+1) gt bri then table(1,counter)=mass1 else table(1,counter)=1.0
    	time=tra(timecol, wp(i):wp(i+1)-1)
    	xavg=float(mean(double(tra(xcol,wp(i):wp(i+1)-1))))
	;temp(counter, round(time))=round((tra(xcol,wp(i):wp(i+1)-1)-xavg)/res)*res
    	temp(counter, round(time))=(tra(xcol,wp(i):wp(i+1)-1)-xavg)
	table(0,counter)=tra(idcol,wp(i)) ;track ID
	table(2,counter)=xavg
	counter=counter+1
	endif
endfor

for i=0, nwp-1 do begin;npart-1 do begin change by PJY
	life=wp(i+1)-wp(i)-1
	if life ge lifetime then begin
        if tra(bcol, wp(i)+1) gt bri then table(1,counter)=mass1 else table(1,counter)=1.0
	time=tra(timecol, wp(i):wp(i+1)-1)
	yavg=float(mean(double(tra(ycol,wp(i):wp(i+1)-1))))
	;temp(counter, round(time))=round((tra(ycol,wp(i):wp(i+1)-1)-yavg)/res)*res
	temp(counter, round(time))=(tra(ycol,wp(i):wp(i+1)-1)-yavg)
	table(0,counter)=tra(idcol,wp(i)) ;track ID    
	table(2,counter)=yavg
    	counter=counter+1
	endif
endfor

write_gdf, table(*,0:counter-1), 'X0XnY0Yn.gdf'

dmx=temp(0:counter-1,*)  ;; displacment matrix
write_gdf, dmx, 'disp_X0XnY0Yn.gdf'
dcmx=mgcorrelate(dmx,/covariance)      ;; construct displacement correlation matrix-memory and cpu intensive
write_gdf, dcmx, 'disp_corr_X0XnY0Yn.gdf'
inv=invert(dcmx)
;dcmx=0
write_gdf, inv, 'km.gdf'
inv=(inv+transpose(inv))/2.0

m=fltarr(counter)
m(*)=table(1,0:counter-1) ;set up masses of particles 
mmask=m#m
mmask=sqrt(mmask)
inv=inv/mmask
inv=(inv+transpose(inv))/2.0

write_gdf,inv,'km_mi.gdf'

;if not keyword_set(ev) then eigen=eigenql(float(inv)) else begin
  eigen=eigenql(float(inv), eigenvectors=egv)  ;; will increases memory usage
  write_gdf, egv, 'eigenvectors.gdf'      ;; the same size as the matrix
  write_gdf,eigen,'eigenvalues.gdf'
  ;write_gdf,fucked,'failed_eigenvalues.gdf'
;endelse

return, dmx
end

