clear
tic

%%%%%%%%%%%%%%%%%%
%%% PARAMETERS %%%
%%%%%%%%%%%%%%%%%%
beta = 0.96; %discount factor
alpha = 1/3; %capital share
delta = 0.025; %depreciation rate 
eta = 1; %inverse Frisch elasticity
xi = 1; %labor supply parameter
ro = 0.95; %TFP autocorrelation

%%%%%%%%%%
%%% SS %%%
%%%%%%%%%%
x = fsolve(@(x) [x(1)-x(2)^alpha*x(3)^(1-alpha);
x(4)-x(1)+delta*x(2);
xi/eta*x(3)^(eta-1)*x(4)-(1-alpha)*x(1)/x(3);
1-beta*(1-delta+alpha*x(1)/x(2))],ones(4,1)); 
y_SS = x(1);
k_SS = x(2);
l_SS = x(3);
c_SS = x(4);


%%%%%%%%%%%%%%%%%%%%%%
%%% BLANCHARD-KAHN %%%
%%%%%%%%%%%%%%%%%%%%%%
phi = alpha*y_SS/k_SS/(1-delta+alpha*y_SS/k_SS); %to simplify notation
A = [0 0 0 0;0 1 0 0; 0 0 1 0; -phi*(1-alpha) -phi phi*(1-alpha) 1];
B = [1-alpha-eta 1 alpha -1; 0 ro 0 0; (1-alpha)*y_SS/k_SS y_SS/k_SS 1-delta+alpha*y_SS/k_SS -c_SS/k_SS;0 0 0 1];
AA = zeros(4,4);
BB = zeros(4,4);
for i = 2:4
    prop1 = B(i,1)/B(1,1);
    prop2 = A(i,1)/B(1,1);
    BB(i,:) = B(i,:)-prop1*B(1,:);
    AA(i,:) = A(i,:)-prop2*B(1,:);
end
AA = AA(2:4,2:4);
BB = BB(2:4,2:4);
D = AA\BB;
[V,L] = eig(D);
Q = inv(V);
j = find(diag(L) > 1); %number of eigenvalue > 1
coint = Q(j,:);
DD = zeros(2,3);
for i = 1:2
    prop3 = D(i,3)/coint(1,3);
    DD(i,:) = D(i,:)-prop3*coint; 
end
P = [DD(:,1:2);-coint(1,1:2)/coint(1,3)];
P(4,:) = -B(1,2:3)/B(1,1)-B(1,4)/B(1,1)*P(3,:);


%%%%%%%%%%%
%%% IRF %%%
%%%%%%%%%%%
T = 100; %IRF horizon
z0 = 1;
state(:,1) = [z0;0];
control = zeros(2,T);
output = zeros(1,T);
invest = zeros(1,T);
for t = 1:T
    control(:,t) = P(3:4,:)*state(:,t); %Consumption and labor
    state(:,t+1) = P(1:2,:)*state(:,t); %Future capital and TFP 
    output(1,t) = state(1,t)+alpha*state(2,t)+(1-alpha)*control(2,t);
    invest(1,t) = state(2,t+1)-(1-delta)*state(2,t); 
end
h = (1:T);


subplot(3,2,1)
plot(h,state(1,1:T))
title('TFP')

subplot(3,2,2)
plot(h,state(2,1:T))
title('Capital')

subplot(3,2,3)
plot(h,control(1,1:T))
title('Consumption')

subplot(3,2,4)
plot(h,control(2,1:T))
title('Labor')

subplot(3,2,5)
plot(h,output(1,1:T))
title('Output')

subplot(3,2,6)
plot(h,invest(1,1:T))
title('Investment')

toc