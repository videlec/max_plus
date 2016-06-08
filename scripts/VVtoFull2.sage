from max_plus import *

d=2#dim
N=(lcm(range(1,d+1)))

a,b=symbolic_max_plus_matrices(d,2)
D = {'0': a^2,'1': b^2}


Id=[]
for l in range(11,12):
	f=open('max_plus/vv_identities/'+str(d)+'_'+str(l))
	text=f.read()
	f.close()
	lines=text.split('\n')
	lines=[l for l in lines if l]
	Id=[tuple(l.split(' ')) for l in lines]+Id

Idf=[]
Idpf=[]
for u,v in Id:
	if (prod(D[i] for i in u))==(prod(D[i] for i in v)):
		Idf=Idf+[(u,v)]
	else:
		Idpf=Idpf+[(u,v)]


