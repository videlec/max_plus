from max_plus import *
from max_plus.combinat import runs

#sage: f=open("tmp_full_2_23","w")
#sage: for I in l:
#....:     f.write(''.join(map(str,I[0])))
#....:     f.write(' ')
#....:     f.write(''.join(map(str,I[1])))
#....:     f.write('\n')
#....:     f.flush()
#....:     
#sage: f.close()


f=open('/home/glenn/Dropbox/Bureau/Glenn/pdfMaths/MPSemiGpes/max_plus/vv_identities/3_23')
N=1000000 #Nb tests
d=3
K=100
p=0.01

for k in range(17,24):
	f=open("2_"+str(k))
	text=f.read()
	f.close()
	lines=text.split('\n')# separe aux saut ligne (les efface) et fait une liste des intervalles
	lines=[l for l in lines if l]
	Id[k]=[tuple(l.split(' ')) for l in lines if l]
	Idc[k]=[(u.replace('0','2').replace('1','0').replace('2','1'),v.replace('0','2').replace('1','0').replace('2','1')) for (v,u) in Id[k]]+Id[k]
	aId[k]=set([('0'+u,'0'+v) for u,v in Idc[k]])
	Ida[k]=set([(u+'0',v+'0') for u,v in Id[k]]+[(u+'1',v+'1') for u,v in Id[k]])

for k in range(18,24):
	idm=set(Id[k]).difference(aId[k-1].union(Ida[k-1]))
	f=open('2_'+str(k)+'m',"w")
	for (u,v) in idm:
		f.write(u)
		f.write(' ')
		f.write(v)
		f.write('\n')
	f.close()

for k in range(18,24):
	idm=(aId[k-1].union(Ida[k-1])).difference(set(Id[k]))
	print(idm)








    
    
rC=[]

for l in lines:
	u,v=l.split(' ')
	n=len(u)
	T=True
	i=1
	while i<N and T:
		a=random_integer_max_plus_matrix(d,0,K,p)
		b=random_integer_max_plus_matrix(d,0,K,p)
		a6=a**6; b6=b**6
		m=[a6,b6]
		T=prod(m[int(u[k])] for k in range(n))==prod(m[int(v[k])] for k in range(n))
		i=i+1
	if T:
		rC=rC+[[u,v]]


