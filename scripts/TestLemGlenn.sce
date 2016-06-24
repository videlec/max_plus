d=5;//taille matrices
st='';
for i=1:d
  st=strcat([st,'%i,']);
end
N=factorial(d);//Puissance

//n=10^7;//Nb tests
//K=10^4;//Size of entries

//Mot
u=[0,1,0,1,0,1,1,0,1,0,1,0,0,1,1,0,0,1,0,1,0,1,0,1,0,1,0,1,1,0,1,0,1,0,0,1,1,0,0,1,0,1,1,0,0,1,0,1,0,1,1,0,1,0,1,0,0,1,1,0,0,1,0,1];
l=size(u,'*');
L=zeros(1,l);
Le=zeros(1,l);
Lee=zeros(1,l);

// Let's go
fidw=mopen('wrank.txt','w');
fidt=mopen('addtr.txt','w');
mfprintf(fidt,'u=')
for i=1:size(u,'*')
	mfprintf(fidt,'%i',u(i))
end
mfprintf(fidt,'\n')

for j=1:n 

A=round(K*rand(d,d));
B0=(-K/2)*#(A);
A=round(K*rand(d,d));
A0=(-K/2)*#(A);
A=A0^N;
B=B0^N;

P=%eye(d,d);
b=%T;
k=0;
D=%ones(d,1);

while k<l & b
    k=k+1;
    if u(k)==0
        C=A;
    else
        C=B;
    end
    P=P*C;
    for i=1:d
      D(i)=D(i)*C(i,i);
    end
    b=and(D==diag(P));
end

L(k)=L(k)+1;

if k==l&(or(A(1,1)*%ones(d,1)<>diag(A))|or(B(1,1)*%ones(d,1)<>diag(B)))
   mfprintf(fidt,st,plustimes(A0))
   mfprintf(fidt,';')
   mfprintf(fidt,st,plustimes(B0))
   mfprintf(fidt,'\n')
end

if (k<l)&size(weakbasis(P),'*')==d*d
    A1=A0;
    B1=B0;
    Le(k)=Le(k)+1;
    if size(weakbasis(P'),'*')==d*d
    	Lee(k)=Lee(k)+1;
        mfprintf(fidw,st,plustimes(P))
        mfprintf(fidw,'\n')
    end

end

end
//L
mclose(fidw);
mclose(fidt);
