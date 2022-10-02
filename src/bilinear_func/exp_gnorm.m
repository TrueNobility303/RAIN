rng(42);

% f(x,y) = xy
% comparing RAIN, EG, EAG, R-EG

% tuned examples 
% 1. n = 900, tau = 0.1, sigma = 0.001
% 2. n =  1600, tau = 0.05, sigma = 0.005
% 3. n = 4000, tau = 0.01, sigma = 0.01
% 4. n = 6000, tau = 0.005, sigma = 0.02

n = 900;
tau = 0.1;
z0 = [10,10]';
sigma = 0.001;

REG_flow = zeros(n,2);
REG_gnorm = zeros(n,1);
z = z0;
lambda = 0.1;
for i = 1:n
    w = z - tau * [z(2),-z(1)]' + tau * lambda * (z0 - z) + normrnd(0,sigma,[2,1]);
    z = z - tau * [w(2),-w(1)]' + tau * lambda * (z0 - w) + normrnd(0,sigma,[2,1]);
    REG_flow(i,:) = z;
    REG_gnorm(i) = norm(z);
end

EAG_flow = zeros(n,2);
EAG_gnorm = zeros(n,1);
z = z0;
for i = 1:n
    w = z - tau * [z(2),-z(1)]' + 1 / (i+2) * (z0 - z) + normrnd(0,sigma,[2,1]);
    z = z - tau * [w(2),-w(1)]' + 1 / (i+2) * (z0 - z) + normrnd(0,sigma,[2,1]);
    EAG_flow(i,:) = z;
    EAG_gnorm(i) = norm(z);
end

RAIN_flow = zeros(n,2);
RAIN_gnorm = zeros(n,1);
z = z0;
gamma = 0.001; 
for i = 1:n
    w = z - tau * [z(2),-z(1)]' + normrnd(0,sigma,[2,1]);
    for  j = 1:i
        w = w + tau * lambda * gamma * (1 + gamma)^j * (RAIN_flow(j,:)' - z);
    end
    z = z - tau * [w(2),-w(1)]' + normrnd(0,sigma,[2,1]); 
    for  j = 1:i
        z = z + tau * lambda  * gamma * (1+ gamma)^j * (RAIN_flow(j,:)' - w) ;
    end
    RAIN_flow(i,:) = z;
    RAIN_gnorm(i) = norm(z);
end

EG_flow = zeros(n,2);
EG_gnorm = zeros(n,1);
z = z0;
for i = 1:n
    w = z - tau * [z(2),-z(1)]' + normrnd(0,sigma,[2,1]);
    z = z - tau * [w(2),-w(1)]' + normrnd(0,sigma,[2,1]);
    EG_flow(i,:) = z;
    EG_gnorm(i) = norm(z);
end

figure(2);
semilogy(EG_gnorm(1:n),'-.','linewidth',4);
hold on;
semilogy(REG_gnorm(1:n),'--','linewidth',4);
hold on;
semilogy(EAG_gnorm(1:n),':','linewidth',4);
hold on;
semilogy(RAIN_gnorm(1:n),'-','linewidth',4);
hold on;
legend('SEG','R-SEG','SEAG','RAIN','fontsize',15);
xlabel('iterations','fontsize',20);
ylabel('gradient norm','fontsize',20);