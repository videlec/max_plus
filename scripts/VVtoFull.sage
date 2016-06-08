from max_plus import *

d=3#dim
N=(lcm(range(1,d+1)))

k=1000000 #Nb paires
K=10000#Taille coeffs
p=0.01#Proba -\infty

Id=[]
for l in range(22,24):
	f=open('max_plus/vv_identities/'+str(d)+'_'+str(l))
	text=f.read()
	f.close()
	lines=text.split('\n')
	lines=[l for l in lines if l]
	Id=[tuple(l.split(' ')) for l in lines]+Id


elts = [((random_integer_max_plus_matrix(d,0,K,p))^N, (random_integer_max_plus_matrix(d,0,K,p))^N) for _ in range(k)]


Id2=[(u,v) for u,v in Id if is_relation(tuple([int(i) for i in u]), tuple([int(i) for i in v]), elts, False)]

print(d,K,k,p,Id2)
