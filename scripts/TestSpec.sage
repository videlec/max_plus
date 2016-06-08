from max_plus import *

d=2#dim
N=(lcm(range(1,d+1)))

k=10^8 #Nb paires
K=10000#Taille coeffs
p=0#Proba -\infty (=0 obligatoire)

T=True
l=0
while (l<k)&T:
	l=l+1
	A=(random_integer_max_plus_matrix(d,-K,K,p))
	B=(random_integer_max_plus_matrix(d,-K,K,p))
	A2=A*A;B2=B*B;
	X=A2*B2;Y=B2*A2;
	P=X*Y
	S=Y*X
	U=P*X*S
	V=P*Y*S
	Ul=U.list();Vl=V.list()
	u1=Ul[0];v1=Vl[0]
	T=[u1+z for z in Vl]==[v1+z for z in Ul]

print(l,A,B,U,V)
