#Vector creation
#Vector-string
x=[1,2,3]
x=[1 2 3]
#Vector-column
x=[1;2;3]
#Sequence-vector
x=(1:1:10)
x=(1:3:10)
#Matrix creation
x=eye(5)
x=ones(3,2)
x=zeros(4,3)
x=[1 2;3 4]
xx=repmat(x,3,1)
xx=kron([1;1;1],x)
R=reshape(xx,3,4)
R=reshape(xx,3,[])
#Changing elements in matrix
A=zeros(3,3)
B=eye(2)
A(1:2,2:3)=B
#Matrix (vector) operations
B=ones(3,3)
A+B
A*B
A^2
xx=inv(x)
xx=x\eye(2)
#By-element operations
A.*B
A.^2
log(B)
exp(B)

#Simple code
clear;
tic;
S=[1 2;3 4];
SS=S\eye(2);
M=eye(2,3);
toc;

#if-statements

x=100;
if (x>100 && x<1000) || x==1001
  y=1;
elseif x==100
  y=0;
else
  y=-1;
end
y

#loops
n = 50;
A = zeros(n,n);
for i = 1:n
  A(i,i)=1;
end
A

x = 0;
for i =1:n
  x=x+i^2;
end
x

x = 0;
i = 1;
while i<51
  x=x+i^2;
  i=i+1;
end
x

#Some other usefull things
doc eye;

Z=inf;

A = [1 2 3; 4 5 6];
Y = max(A,3)
Y = max(A, [], 1)
[Y, I] = max(A, [], 1)

W=[3 7 8 11 2 9 3 5 4 1];
Q=find(W>3)

#Interpolation
X=[0.5,1,2,3];
Y=[1,1.5,3,0];
Z=[1.5,1.6,2.2,2.9];
Yz=interp1(X,Y,Z, 'linear')
Yz=interp1(X,Y,Z, 'cubic')
Yz=interp1(X,Y,Z, 'spline')

#Solving equations
fsolve(@(x) x^2-5*x+6, 2.5)
fsolve(@(x) [x^2-5*x+6;x-2], 2.5)

#-2x^2 + 3xy + 4sin(y) = 6
#3x^2 - 2xy^2 + 3cos(x) = -4

fsolve(@(x) [-2*x(1)^2 + 3*x(1)*x(2) + 4*sin(x(2)) - 6;
3*x(1)^2 - 2*x(1)*x(2)^2 + 3*cos(x(1)) + 4], ones(1,2))

#Optimizing
[x, fval, info]=fminsearch(@(x) 3*x(1)^2+5*(x(2)-5)^2+2*(x(3)+2)^2, ones(1,3))
[x, fval, info]=fminsearch(@(x) 3*x(1)^2+5*(x(2)-5)^2+2*(x(3)+2)^3, ones(1,3))
[x, fval, info]=fminsearch(@(x) sin(x), 3)