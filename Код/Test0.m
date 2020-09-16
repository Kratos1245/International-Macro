#Interpolation
X=[0.5,1,2,3];
Y=[1,1.5,3,0];
Z = (0.5:0.25:3)
Yz = interp1(X, Y, Z, 'cubic')
plot(Z, Yz)


