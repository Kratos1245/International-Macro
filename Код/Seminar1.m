clear
tic

%%%%%%%%%%%%%%%%%%
%%% PARAMETERS %%%
%%%%%%%%%%%%%%%%%%
beta = 0.96; %discount factor
alpha = 1/3; %elasticity of production wrt capital
A = 1; %TFP
sigma = 1; %1/elasticity of intertemporal substitution
delta = 0.1; %depreciation rate

%%%%%%%%%%
%%% SS %%%
%%%%%%%%%%
k_SS = ((1/beta-1+delta)/(alpha*A))^(1/(alpha-1)); %capital in SS
c_SS = A*k_SS^alpha-delta*k_SS; %consumption in SS
c_min = 0
c_max = c_SS

%%%%%%%%%%%%%%%%%%%
%%% SIMULATIONS %%%
%%%%%%%%%%%%%%%%%%%
L = 100; %number of periods of simulation
k(1) = 1;
c(1) = 0.603; %0.1, 0.5, 0.7, 0.6, 0.603

%Simple iteration (do not run together with the iteration with conditions)

for t = 2:L
    k(t) = A*k(t-1)^alpha + (1-delta)*k(t-1)-c(t-1)
    c(t) = beta*(1-delta+alpha*A*k(t)^(alpha-1))*c(t-1)
end

%Iteration with conditions (do not run together with the simple iteration)
t=2
while t < L
  if abs(c_SS-c(t-1))>0.001 || abs(k_SS-k(t-1))>0.001
    k(t) = A*k(t-1)^alpha + (1-delta)*k(t-1)-c(t-1)
    c(t) = beta*(1-delta+alpha*A*k(t)^(alpha-1))*c(t-1)
    if k(t)-k(t-1)<0
      c_max = c(1)
      c(1)=(c(1)+c_min)/2
      t=2
    elseif c(t)-c(t-1)<0
      c_min = c(1)
      c(1)=(c(1)+c_max)/2
      t=2
    else
      t=t+1
    end
  else
    break,
  end
end

%%%%%%%%%%%%%%%
%%% FIGURES %%%
%%%%%%%%%%%%%%%

plot(k,c)
title('Seddle path')
xlabel('Capital')
ylabel('Consumption')
legend('Seddle path','Location','Best')

%%%%%%%%%%%%%%%
%%% SHOCK %%%%%
%%%%%%%%%%%%%%%

A_new = 1.1
k_SS_new = (( 1/beta-1+delta)/(alpha*A_new))^(1/(alpha-1)); %capital in SS
c_SS_new = A_new*k_SS_new^alpha-delta*k_SS_new; %consumption in SS

c_new(1) = c_SS+0.085 %-0.01 +0.01 +0.08
k_new(1) = A_new*k_SS^alpha + (1-delta)*k_SS-c_new(1)

L = 100

for t = 2:L
    k_new(t) = A_new*k_new(t-1)^alpha + (1-delta)*k_new(t-1)-c_new(t-1)
    c_new(t) = beta*(1-delta+alpha*A_new*k_new(t)^(alpha-1))*c_new(t-1)
end

plot(k_new,c_new)
title('Seddle path')
xlabel('Capital')
ylabel('Consumption')
legend('Seddle path','Location','Best')

%%%%%%%%%%%%%%%%%%%%%%
%%% BLANCHARD-KAHN %%%
%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%
%%% SS %%%
%%%%%%%%%%
k_SS = ((1/beta-1+delta)/(alpha*A))^(1/(alpha-1)); %capital in SS
c_SS = A*k_SS^alpha-delta*k_SS; %consumption in SS
M = [beta*alpha*(1-alpha)*A*k_SS^(alpha-1) sigma; 1 0];
B = [0 sigma; 1/beta -c_SS/k_SS];
D = M\B;
[V,L] = eig(D); 
Q = inv(V);

toc