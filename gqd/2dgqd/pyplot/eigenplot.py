#import moduls
import numpy as np
import matplotlib.pyplot as plt
import sys, linecache

num = int(sys.argv[1])
eigenvalue = "/Users/guqiangqiang/github/Eigenproblem/gqd/2dgqd/eigva.dat"
eig = linecache.getline(eigenvalue,num)
print eig

eigenvector = "/Users/guqiangqiang/github/Eigenproblem/gqd/2dgqd/eigve.dat"
eiv = linecache.getline(eigenvector,num)
eiv = eiv.split()
M = 50
N = 50
W = 8
L = 8
Y = eiv 
k = 0
fig = plt.figure(figsize=(8,6.5))
ZZ = np.zeros((N-1,M-1))
for j in range (0,M-1):
    for i in range (0,N-1):
        ZZ[i][j]=float (Y[k])
        k=k+1
y=np.linspace(-W/2,W/2,N-1)
x=np.linspace(-L/2,L/2,M-1)
[X,Y]=np.meshgrid(x,y)

plt.contourf(X, Y, ZZ, 8, alpha=0.65, cmap=plt.cm.jet)
#C = plt.contour(X, Y, ZZ, 8, colors='black', linewidth=.5)
#plt.clabel(C, inline=1, fontsize=10)
plt.colorbar()
plt.show()