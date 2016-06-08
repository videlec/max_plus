f=open("max_plus/vv_identities/3_22")
text=f.read()
f.close()
lines=text.split('\n')
lines=[l for l in lines if l]
IdVV322=[tuple(l.split(' ')) for l in lines if l]
IdVV322c=[(u.replace('0','2').replace('1','0').replace('2','1'),v.replace('0','2').replace('1','0').replace('2','1')) for (v,u) in IdVV322]+IdVV322
IdVV322c=[(u[::-1],v[::-1]) for u,v in IdVV322c]+IdVV322c
IdVV322c=[(v,u) for u,v in IdVV322c]+IdVV322c


f=open("max_plus/sv_identities/3_11")
text=f.read()
f.close()
lines=text.split('\n')
lines=[l for l in lines if l]
IdSV311=[tuple(l.split(' ')) for l in lines if l]
IdSV311c=[(u.replace('0','2').replace('1','0').replace('2','1'),v.replace('0','2').replace('1','0').replace('2','1')) for (v,u) in IdSV311]+IdSV311
IdSV311c=[(u[::-1],v[::-1]) for u,v in IdSV311c]+IdSV311c
IdSV311c=[(v,u) for u,v in IdSV311c]+IdSV311c

IdSub311c=[(u.replace('1','2').replace('0','01').replace('2','10'),v.replace('1','2').replace('0','01').replace('2','10')) for u,v in IdSV311c]



print(len(set(IdSub311c).difference(set(IdVV322c))))

print(len(set(IdVV322c).difference(set(IdSub311c))))



