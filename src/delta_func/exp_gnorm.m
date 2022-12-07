

% a wost case "delta,nu" function in 
% paper:Performance of first-order methods for smooth convex minimization: A novel approach

rng(42);
nu = 5e-5;
delta = 1e-2;

d = 100;
tau = 1;
n = 6000;
z0 = randn(2*d,1);
sigma = 0.02;

REG_flow = zeros(n,2*d);
REG_gnorm = zeros(n,1);
z = z0;
lambda = 0.01;
for i = 1:n
    gz = full_grad(nu,delta,z);
    w = z - tau * gz + tau * lambda * (z0 - z) + normrnd(0,sigma,[2*d,1]);
    gw = full_grad(nu,delta,w);
    z = z - tau * gw+ tau * lambda * (z0 - w) + normrnd(0,sigma,[2*d,1]);
    REG_flow(i,:) = z';
    REG_gnorm(i) = norm(gz);
end

EAG_flow = zeros(n,2*d);
EAG_gnorm = zeros(n,1);
z = z0;
for i = 1:n
    gz = full_grad(nu,delta,z);
    w = z - tau * gz + 1 / (i+2) * (z0 - z) + normrnd(0,sigma,[2*d,1]);
    gw = full_grad(nu,delta,w);
    z = z - tau * gw + 1 / (i+2) * (z0 - z) + normrnd(0,sigma,[2*d,1]);
    EAG_flow(i,:) = z';
    EAG_gnorm(i) = norm(gz);
end

RAIN_flow = zeros(n,2*d);
RAIN_gnorm = zeros(n,1);
z = z0;
gamma = 0.0001; 
for i = 1:n
    gz = full_grad(nu,delta,z);
    w = z - tau * gz+ normrnd(0,sigma,[2*d,1]);
    for  j = 1:i
        w = w + tau * lambda * gamma * (1 + gamma)^j * (RAIN_flow(j,:)' - z);
    end
    gw = full_grad(nu,delta,w);
    z = z - tau * gw + normrnd(0,sigma,[2*d,1]); 
    for  j = 1:i
        z = z + tau * lambda  * gamma * (1+ gamma)^j * (RAIN_flow(j,:)' - w) ;
    end
    RAIN_flow(i,:) = z';
    RAIN_gnorm(i) = norm(gz);
end

EG_flow = zeros(n,2*d);
EG_gnorm = zeros(n,1);
z = z0;
for i = 1:n
    gz = full_grad(nu,delta,z);
    w = z - tau * gz + normrnd(0,sigma,[2*d,1]);
    gw = full_grad(nu,delta,w);
    z = z - tau * gw + normrnd(0,sigma,[2*d,1]);
    EG_flow(i,:) = z';
    EG_gnorm(i) = norm(gz);
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

function g = full_grad(eps,delta,z)
    d = length(z) / 2;
    x = z(1:d);
    y = z(d+1:end);
    gx = (1-delta) * sub_grad(eps,x) + delta * y;
    gy = (1-delta) * sub_grad(eps,y) - delta * x;
    g = [gx;gy];
end

function g = sub_grad(eps,u)
    g = u;
    g(u>eps) = eps;
    g(u<eps) = -eps;
end
