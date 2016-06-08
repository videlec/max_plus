from max_plus import *

d=3#dim
N=(lcm(range(1,d+1)))*2*(d-1)

k=10^5 #Nb paires
K=10^3#Taille coeffs
p=0#Proba -\infty (=0 obligatoire)

elts = [((random_integer_max_plus_matrix(d,0,K,p))^N, (random_integer_max_plus_matrix(d,0,K,p))^N) for _ in range(k)]




Id=[]
for l in range(22,23):
	f=open('max_plus/vv_identities/'+str(d)+'_'+str(l))
	text=f.read()
	f.close()
	lines=text.split('\n')
	lines=[l for l in lines if l]
	Id=[tuple(l.split(' ')) for l in lines]+Id

Idf=[]
Idpf=[]
Idt=[]
for u,v in Id:
	FI=True;PI=True
	for a,b in elts:
		D= {'0': a,'1': b}
		U=(prod(D[i] for i in u)); V=(prod(D[i] for i in v))
		if U<>V:
			FI=False;
			Ul=U.list();Vl=V.list()
			x=Ul[0];y=Vl[0]
			if [ul+y for ul in Ul]<>[ul+x for ul in Vl]:
				PI=False
	if FI:
		Idf=Idf+[(u,v)]
	else:
		if PI:
			Idpf=Idpf+[(u,v)]
		else:
			Idt=Idt+[(u,v)]



