#!/usr/bin/env/python

import numpy as np
import math as m
import matplotlib.pyplot as plt


# extract high symmetry K from KPOINTS
kp = open('KPOINTS', 'r').readlines()
k1 = filter(lambda num: num in '1234567890. ', kp[4]).split()
k2 = filter(lambda num: num in '1234567890. ', kp[5]).split()
k3 = filter(lambda num: num in '1234567890. ', kp[8]).split()
kp1 = [ float(i) for i in k1 ]
kp2 = [ float(i) for i in k2 ]
kp3 = [ float(i) for i in k3 ]

num = int(kp[1])

# calculate the distance between high symmetry points
d1 = m.sqrt((kp2[0]-kp1[0])**2 + (kp2[1]-kp1[1])**2 + (kp2[2]-kp1[2])**2)
d2 = m.sqrt((kp3[0]-kp2[0])**2 + (kp3[1]-kp2[1])**2 + (kp3[2]-kp2[2])**2)
d3 = m.sqrt((kp3[0]-kp1[0])**2 + (kp3[1]-kp1[1])**2 + (kp3[2]-kp1[2])**2)

d = d1 + d2 + d3
sd1 = d1/d
sd2 = d2/d
sd3 = d3/d

s1 = sd1/num
s2 = sd2/num
s3 = sd3/num

xl = []
for i in range(num):
    xl.append(i*s1)

for i in range(num):
    xl.append(sd1+i*s2)

for i in range(num):
    xl.append(sd1+sd2+i*s3)

#x.append(1)
x = np.array(xl)

#kp.close()

# pick up KPONTS band energetics ~ projection parameter
pro = open('PROCAR','r').readlines()

# get total number of k-points, bands and ions:
info = filter(lambda num: num in '1234567890. ', pro[1]).split()

k = []
# screen k coordinates and band energies and its ratio for MoS2
for i in pro:
    if 'k-point ' in i:
        for j in filter(lambda num: num in '1234567890. ',i[19:51]).split():
            k.append(float(j))
            print j
    elif 'band ' in i:
        k.append(float(i[18:30]))
        print(float(i[18:30]))
        mos = float(pro[pro.index(i)+3][-6:-1]) + float(pro[pro.index(i)+4][-6:-1]) + float(pro[pro.index(i)+5][-6:-1])
        print mos
        wse = float(pro[pro.index(i)+6][-6:-1]) + float(pro[pro.index(i)+7][-6:-1]) + float(pro[pro.index(i)+8][-6:-1])
        print wse
        r = mos/(mos + wse)
        print r
        k.append(r)
#        k.append((float(pro[pro.index(i)+3][-6:-1]) + float(pro[pro.index(i)+4][-6:-1]) + float(pro[pro.index(i)+5][-6:-1]))/(float(pro[pro.index(i)+3][-6:-1]) + float(pro[pro.index(i)+4][-6:-1]) + float(pro[pro.index(i)+5][-6:-1]) + float(pro[pro.index(i)+6][-6:-1]) + float(pro[pro.index(i)+7][-6:-1]) + float(pro[pro.index(i)+8][-6:-1])))

#print(k)
# list to numpy array
d = np.array(k)
# reshape data
data = d.reshape(int(info[0]), int(info[1])*2+3)

# write data
np.savetxt('data', data)


fig, ax = plt.subplots()
xmin, xmax = x.min(), x.max

# plot high symmetric k.
for i in [d1, d1+d2]:
    ax.axvline(x=i, color='k', lw=1.0, ls=':')

# plot scatter
for i in range(int(info[1])):
    plt.scatter(x, k[:,2*i+4], c=k[:,2*i+5], vmin=xmin, vmax=xmax, s=20, cmap=cm)

plt.colorbar()

ax.axhline(y=0, color='r', ls='--', lw=1.0)

ax.set_xlim((xmin, xmax))
ax.set_ylim((-2,3))

pos = [0,] + [d1,] + [d1+d2,] + [1,]
ax.set_xticks(pos)
print pos

ax.set_xticklabels([r'$\Gamma$','M','K',r'$\Gamma$'], pos,fontsize=15)
ax.set_ylabel('Energy [eV]', fontsize=14)

for sp in ax.spines:
    ax.spines[sp].set_linewidth(1.5)

plt.savefig('band.png', dpi=360)
plt.show()
