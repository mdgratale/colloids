FUNCTION mgfourierdisplacementfix,dis,freq
;;dis=displacement matrix,freq=frequency to cut out

np=n_elements(dis[*,0])
nt=n_elements(dis[0,*])

f=fft(dis,dimension=2)
q=findgen(nt)/nt

w=where(q gt 2e-3 and q lt 1-2e-3,nw,complement=z)
h=f[*,w]
k=q[w]

v=where(k mod freq gt 0.05*freq and k mod freq lt 0.95*freq,nv,complement=b)
h[*,b]=0.0

newf=[[f[*,z]],[h]]
newq=[q[z],k]

finalf=newf[*,sort(newq)]

newdis=fft(finalf,dimension=2,/inverse)
stop
return,newdis
END