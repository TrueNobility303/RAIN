rng(42);


% Hardness function from "Fast Extra Gradient Methods for Smooth Structured Nonconvex-Nonconcave Minimax Problems" , NeurIPS 21
% f(x,y) = -1/6 * x^2 + 2 \sqrt{2} / 3 xy + 1/6 * y^2.
% For this function, we have rho = -1/3, L = 1

% comparing RAIN++ and SEG+  

n = 1000;
z0 = [10,10]';
sigma = 0.005;

rho = - 1 / 3;
L = 1;
SEG_SFO = zeros(n,1);
SEG_gnorm = zeros(n,1);
z = z0;
bar_w = z0;

tau = 1 / (4 * sqrt(2) * L);

for i = 1:n
    w = z - tau * (grad(z, rho, L) + normrnd(0,sigma,[2,1]));
    z = z - tau / 2 * (grad(w, rho, L) + normrnd(0,sigma,[2,1])); 
    SEG_SFO(i) = 2 * i;
    SEG_gnorm(i) = norm(grad(w, rho, L));
end


SFEG_SFO = zeros(n,1);
SFEG_gnorm = zeros(n,1);
z = z0;
for i = 1:n
    gz = grad(z, rho, L) + normrnd(0,sigma,[2,1]);
    w = z  + 1 / i * (z0 - z) - ( 1- 1 / i)  * (1 / L  + 2 * rho) * gz;
    z = z  + 1 / i * (z0 - z) -  1 / L * (grad(w, rho, L) + normrnd(0,sigma,[2,1])) - ( 1- 1 / i) * 2 * rho * gz;
    SFEG_SFO(i) = 2*i;
    SFEG_gnorm(i) = norm(grad(w, rho, L));
end

RAIN_SFO = zeros(n,1);
RAIN_flow = zeros(n,2);
RAIN_gnorm = zeros(n,1);
z = z0;
SFO = 0;
w22 = z0;
gamma = 0.0001;
alpha = 0.5;
tau = 0.5;
for i = 1:n
    [w1,w2] = prox2(z, w22,sigma, alpha, rho, L);
    SFO = SFO + 1;
    J = geornd(0.5);
    if J > 0
        SFO = SFO + 1;
        w12 = 2 * w2 - w1;
    else
        w12 = w1;
    end
    w = z + tau * 2 * (w12- z);
    for  j = 1:i
        w = w + tau * lambda * gamma * (1 + gamma)^j * (RAIN_flow(j,:)' - z);
    end

    [w1,w2] = prox2(w, w12,sigma, alpha, rho, L);
    SFO = SFO + 1;
    J = geornd(0.5);
    if J > 0
        SFO = SFO + 1;
        w22 = 2 * w2 - w1;
    else
        w22 = w1;
    end

    z = z + tau * 2 * (w22 - w);
    for  j = 1:i
        z = z + tau * lambda  * gamma * (1+ gamma)^j * (RAIN_flow(j,:)' - w) ;
    end
    RAIN_flow(i,:) = z;
    RAIN_gnorm(i) = norm(grad(w, rho, L));
    RAIN_SFO(i) = SFO;
end

figure(2);
marker = 1: (n/40) : n;
semilogy(SEG_SFO(1:n), SEG_gnorm(1:n),'--','linewidth',4);
hold on;
semilogy(SFEG_SFO(1:n), SFEG_gnorm(1:n),'-o','linewidth',4, 'Color', 	"#EDB120", 'MarkerIndices', marker, 'MarkerSize', 12);
hold on;
semilogy(RAIN_SFO(1:n), RAIN_gnorm(1:n),'-','linewidth',4, 'Color', "#77AC30");
hold on;
legend('SEG+', 'SFEG', 'RAIN++', 'fontsize',15);
xlabel('#SFO','fontsize',20);
xlim([0,500]);
ylabel('Gradient Norm','fontsize',20);
ax = gca; 
ax.FontSize = 15;

% gradient of the function
function g = grad(z, rho, L)
    gx =  rho * L * L * z(1) + L * sqrt(1 - rho * rho * L * L) * z(2);
    gy =  rho * L * L * z(2) - L * sqrt(1 - rho * rho * L * L) * z(1);
    g = [gx;gy];
end

% k step of SEG on the proximal sub-problem
function w = prox(z, w, k, sigma, alpha, rho, L)
    for i=1:k
        w12 = w -  alpha * (grad(w, rho, L) + 2 * (w - z) + normrnd(0,sigma,[2,1])); 
        w = w - alpha * (grad(w12, rho, L) + 2 * (w12 - z) + normrnd(0,sigma,[2,1]));
    end
end

% 2 step of SEG on the proximal sub-problem
function [w1,w2] = prox2(z, w, sigma, alpha, rho, L)
    w1 = w -  alpha * (grad(w, rho, L) + 2 * (w - z) + normrnd(0,sigma,[2,1])); 
    w = w - alpha * (grad(w1, rho, L) + 2 * (w1 - z) + normrnd(0,sigma,[2,1]));
    w2 = w -  alpha * (grad(w, rho, L) + 2 * (w - z) + normrnd(0,sigma,[2,1]));
end