from max_plus import *
from max_plus.max_plus_int import IntegerMatrixProduct, IntegerMaxPlusMatrix

d=3#dim
N=(lcm(range(1,d+1)))

k=7#10^k #Nb paires
K=10^5#Taille coeffs
p=0.01#Proba -\infty

Id2=[]
for l in range(18,20):
	f=open('../../max_plus/full_identities/'+str(d-1)+'_'+str(l))
	text=f.read()
	f.close()
	lines=text.split('\n')
	lines=[l for l in lines if l]
	Id2=[tuple(l.split(' ')) for l in lines]+Id2
#u,v=Id2[0]

Id3=[]
for l in range(22,24):
	f=open('../../max_plus/vv_identities/'+str(d)+'_'+str(l))
	text=f.read()
	f.close()
	lines=text.split('\n')
	lines=[l for l in lines if l]
	Id3=[tuple(l.split(' ')) for l in lines]+Id3
#s,t=Id3[0]

pN=IntegerMatrixProduct((0,)*N)
M=IntegerMaxPlusMatrix(d,d, [0]*(d*d))
eltsAB = [((random_integer_max_plus_matrix(d,0,K,p))^N, (random_integer_max_plus_matrix(d,0,K,p))^N) for _ in range(10^k)]

Id3f=[]
f=open("tmp_full_"+str(d),"w")

Nb=[0]*k

for u,v in Id2:
	Pu=IntegerMatrixProduct(u)
	Pv=IntegerMatrixProduct(v)
	eltsUV=[]
	for a,b in eltsAB:
		U=Pu(a,b)
		V=Pv(a,b)
		if U<>V:
			eltsUV.append((U,V))

	for s,t in Id3:
		l=1
		T=True
		while l<k:
			if is_relation(tuple([int(i) for i in s]), tuple([int(i) for i in t]), eltsUV[0:(10^(l+1))], False):
				Nb[l]=Nb[l]+1
			else:
				l=k
			l=l+1
			if l==k:
				Id3f.append((s,t,u,v))bloquer
				f.write(s+' '+t' '+u+' '+v)
				f.write('\n')
				f.flush()                
                
               
f.close()

print(d,K,k,p)
print(Nb)
