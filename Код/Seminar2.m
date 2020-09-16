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
c_min = 0;
c_max = c_SS;

%%%%%%%%%%%%%%%%%%%
%%% SIMULATIONS %%%
%%%%%%%%%%%%%%%%%%%
L = 1000; %number of periods of simulation
k(1) = 1;
c(1) = 0.603; %0.1, 0.5, 0.7, 0.6, 0.603

t=2
while t < L
  if abs(c_SS-c(t-1))>0.001 || abs(k_SS-k(t-1))>0.001;
    k(t) = A*k(t-1)^alpha + (1-delta)*k(t-1)-c(t-1);
    c(t) = beta*(1-delta+alpha*A*k(t)^(alpha-1))*c(t-1);
    if k(t)-k(t-1)<0
      c_max = c(1);
      c(1)=(c(1)+c_min)/2;
      t=2;
    elseif c(t)-c(t-1)<0
      c_min = c(1);
      c(1)=(c(1)+c_max)/2;
      t=2;
    else
      t=t+1;
    end
  else
    break,
  end
end

%%%%%%%%%%%%%%%%%%%%%%
%%% BLANCHARD-KAHN %%%
%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%
%%% SS %%%
%%%%%%%%%%
M = [beta*alpha*(1-alpha)*A*k_SS^(alpha-1) sigma; 1 0];
B = [0 sigma; 1/beta -c_SS/k_SS];
D = M\B;
[V,L] = eig(D); 
Q = inv(V);
k_bk(1)=(1-k_SS)/k_SS;
c_bk(1)=-Q(1,1)/Q(1,2)*k_bk(1);
k_bk_real(1)=1;
c_bk_real(1)=c_bk(1)*c_SS+c_SS;


%%%%%%%%%%%%%%%%%%%
%%% SIMULATIONS %%%
%%%%%%%%%%%%%%%%%%%
N=100
num = find(diag(L) > 1);
for t = 2:N
  k_bk_real(t)=A*k_bk_real(t-1)^alpha + (1-delta)*k_bk_real(t-1)-c_bk_real(t-1);
  k_bk(t) = (k_bk_real(t)-k_SS)/k_SS;
  c_bk(t) = -Q(num,1)/Q(num,2)*k_bk(t);
  c_bk_real(t)=c_bk(t)*c_SS+c_SS;
end


%%%%%%%%%%%%%%%
%%% FIGURES %%%
%%%%%%%%%%%%%%%
plot(k, c, k_bk_real, c_bk_real)
title('Seddle path')
xlabel('Capital')
ylabel('Consumption')
legend('Seddle path SA', 'Seddle path BK', 'Location','Best')

%%%%%%%%%%%%%%%%%%%%%%%
%%% TRANSITION PATH %%%
%%%%%%%%%%%%%%%%%%%%%%%
eps = 0.0001;
k_t(1) = k_SS-eps;
k_bk = (k_t(1)-k_SS)/k_SS;
c_bk = -Q(num,1)/Q(num,2)*k_bk;
c_t(1) = c_bk*c_SS+c_SS;
T = 100;
for t = 2:T
    c_t(t) = c_t(t-1)*(beta*(1-delta+A*alpha*k_t(t-1)^(alpha-1)))^(-1);
    k_t(t) = fsolve(@(x) x^alpha+(1-delta)*x-c_t(t)-k_t(t-1),k_t(t-1));
end

plot(k, c, k_bk_real, c_bk_real, k_t, c_t)
title('Seddle path')
xlabel('Capital')
ylabel('Consumption')
legend('Seddle path SA', 'Seddle path BK', 'Seddle path BK backward', 'Location','Best')

%%%%%%%%%%%%
%%% GRID %%%
%%%%%%%%%%%%
k_min = 1;
k_max = k_SS;
n = 1000; %number of elements on the grid
g = (0:1/(n-1):1)';
k_g = k_min+(k_max-k_min)*g; %grid of capital

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SOLVING FOR VALUE FUNCTION %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c_g = A*repmat(k_g,1,n).^alpha+(1-delta)*repmat(k_g,1,n)-repmat(k_g',n,1); %consumption of grid
c_g = max(c_g,0); %consumption of grid
U = log(c_g); %utility on grid
V=[];
V(:,1) = repmat(log(c_SS)/(1-beta),n,1);
[V(:,2), I(:,2)] = max(U+beta*repmat(V(:,1)',n,1),[],2);

t=3;
while max(abs(V(:,size(V,2))-V(:,size(V,2)-1)))>0.000001
  [V(:,t), I(:,t)] = max(U+beta*repmat(V(:,t-1)',n,1),[],2);
  t=t+1;
end

opt_k=k_g(I(:,size(I,2)));
opt_V=V(:,size(V,2));
opt_c=A*k_g.^alpha+(1-delta)*k_g-opt_k;

%%%%%%%%%%%%%%%
%%% FIGURES %%%
%%%%%%%%%%%%%%%

plot(k, c, k_bk_real, c_bk_real, k_t, c_t, k_g, opt_c)
title('Seddle path')
xlabel('Capital')
ylabel('Consumption')
legend('Seddle path SA', 'Seddle path BK', 'Seddle path BK backward', 'Seddle path Bellman', 'Location','northwest')

toc