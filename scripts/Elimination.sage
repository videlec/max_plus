from max_plus import *

d=3#dim
N=(lcm(range(1,d+1)))


T=2#Nb de tours
k=7#10^k #Nb paires A,B
K=10^5#Taille coeffs
p=0.01#Proba -\infty



#Recuperation des candidates
Id8=[]
f=open('max_plus/candidates_full/10^8vv3_22-23(full2_18-20)(a6b6)')
text=f.read()
f.close()
lines=text.split('\n')
lines=[l for l in lines if l]
for l in lines:
	Id8.append(tuple(l.split(' ')))


#S=set(Id8).intersection(set(Id3),set(Id1),set(Id2))
S=Id8


#Tri des candidates (mal fait)
Suv=set()
for s,t,u,v in S:
	Suv.add((u,v))
Suv=list(Suv)
n=len(Suv)

Sst=range(n)
for i in range(n):
	Sst[i]=[]

for s,t,u,v in S:
	for i in range(n):
		if (u,v)==Suv[i]:
			(Sst[i]).append((tuple([int(i) for i in s]), tuple([int(i) for i in t])))


#Eliminations
for _ in range(T):
	eltsAB = [((random_integer_max_plus_matrix(d,0,K,p))^N, (random_integer_max_plus_matrix(d,0,K,p))^N) for _ in range(10^k)]
	for i in range(n):
		u,v=Suv[i]
		Pu=IntegerMatrixProduct(u)
		Pv=IntegerMatrixProduct(v)
		eltsUV=[]
		for a,b in eltsAB:
			U=Pu(a,b)
			V=Pv(a,b)
			if U<>V:
				eltsUV.append((U,V))
		Sst[i]=[st for st in Sst[i] if is_relation(st[0],st[1], eltsUV, False)]

#Res
L=[len(Sst[i]) for i in range(n)]
print(sum(L),'vs',len(S))

if len(S)<sum(L):
	f=open('max_plus/candidates_full/10^8vv3_22-23(full2_18-20)(a6b6)','w')
	for i in range(n):
		u,v=Suv[i]
		for s,t in Sst[i]:
			f.write(''.join(map(str, s))+' '+''.join(map(str, t))+' '+u+' '+v)


