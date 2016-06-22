from max_plus import *

d=5;

f=open('wrank.txt')
text=f.read()
f.close()
lines=text.split('\n')
lines=[l for l in lines if l]
lines=[l[:-1].split(',') for l in lines]
assert len(lines[0])==d*d
matrices=[max_plus_int.IntegerMaxPlusMatrix(d,[int(n) for n in l])for l in lines]
examples=[]
for M in matrices:
	if M.barvinok_rank()==d:
		examples.append(M)



