
% f(x,y) = xy
% comparing RAIN, EG, EAG, R-EG 

n = 1000;
tau = 0.1;
z0 = [10,10]';
sigma = 0.001;

REG_flow = zeros(n,2);
REG_gnorm = zeros(n,1);
z = z0;
for i = 1:n
    w = z - tau * [z(2),-z(1)]' + tau / 10 * (z0 - z) + normrnd(0,sigma,[2,1]);
    z = z - tau * [w(2),-w(1)]' + tau / 10 * (z0 - w) + normrnd(0,sigma,[2,1]);
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
for i = 1:n
    w = z - tau * [z(2),-z(1)]' + normrnd(0,sigma,[2,1]);
    for  j = 1:i
        w = w + tau / 100 * 1.0001^j * (RAIN_flow(j,:)' - z);
    end
    z = z - tau * [w(2),-w(1)]' + normrnd(0,sigma,[2,1]); 
    for  j = 1:i
        z = z + tau / 100  * 1.0001^j * (RAIN_flow(j,:)' - w) ;
    end
    RAIN_flow(i,:) = z;
    RAIN_gnorm(i) = norm(z);
end

EG_flow = zeros(n,2);
EG_gnorm = zeros(n,1);
z = z0;
for i = 1:n
    w = z - tau * [z(2),-z(1)]';
    z = z - tau * [w(2),-w(1)]';
    EG_flow(i,:) = z;
    EG_gnorm(i) = norm(z);
end

figure(1);

subplot(2,2,1);
scatter(EG_flow(:,1), EG_flow(:,2),1);
title('SEG');

subplot(2,2,2);
scatter(REG_flow(:,1), REG_flow(:,2),1);
title('R-SEG');

subplot(2,2,3);
scatter(EAG_flow(:,1), EAG_flow(:,2),1);
title('SEAG');

subplot(2,2,4);
scatter(RAIN_flow(:,1), RAIN_flow(:,2),1);
title('RAIN');


