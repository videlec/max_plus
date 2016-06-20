from max_plus import *

d=3;#Dimension
N=lcm(range(1,d+1)) #Order of matrices
K=1000 #range of entries
n=1000 #Number of tests

#Word to test
p='xxyyx'
u=p+'x'+p+'y'+p
u=u.replace('x','01').replace('y','10')

l=len(u)

L=[0 for _ in range(0,l)]

for _ in range(0,n):
	A=(max_plus_int.random_integer_max_plus_matrix(d,-2*K,K,0))^N;
	B=(max_plus_int.random_integer_max_plus_matrix(d,-2*K,K,0))^N;
	M=[A,B]

	P=integer_max_plus_matrix_identity(d)
	D=[0 for i in range(0,d)]

	k=0;T=True
	while T&(k<l):
		C=M[int(u[k])]
		P=P*C
		D=[D[i]+C[i,i] for i in range(0,d)]
		T=all(D[i]==P[i,i] for i in range(0,d))
		k=k+1
	
	L[k-1]=L[k-1]+1
	if k==l:
		print('Trace add:', A,B,u)
	elif P.barvinok_rank()==d:
		print('VSEx:',A,B,u[0:k+1])


