% computation parameter
N = 5000;  % number of timesteps

% option parameters
S0 = 10;      % initial stock price
K = 10;       % strike price
T = 1;       % expiry time
r = 0.02;      % interest rate
sigma = 0.2;  % volatility
mu = 0.05; %drift

% tree parameters
dt = T/N;         % size of timesteps
t = [0:dt:T];
u = exp(r*dt+sigma*sqrt(dt));  % stock price up ratio
d = exp(r*dt-sigma*sqrt(dt));                  % stock price down ratio
%dscnt = exp(-r*dt);       % discount rate
qu = (1-exp(-sigma*sqrt(dt)))/(exp(sigma*sqrt(dt))-exp(-sigma*sqrt(dt))); % probability (up)
qd = 1 - qu;   % probability (down)

exercise_boundary=[];
exercise_boundary(:,1) = AmericanPut(N, S0, K, T, r, 0.2,mu);
figure(1);
plot(t, exercise_boundary);
ylim([6 10]);
xlabel('time t')
ylabel('S_t')
leg = legend('\sigma=0.2');
title("Exercise Boundary for \sigma=0.2, r=0.02")

exercise_boundary=[];
exercise_boundary(:,1) = AmericanPut(N, S0, K, T, r, 0.1,mu);
exercise_boundary(:,2) = AmericanPut(N, S0, K, T, r, 0.15,mu);
exercise_boundary(:,3) = AmericanPut(N, S0, K, T, r, 0.2,mu);
exercise_boundary(:,4) = AmericanPut(N, S0, K, T, r, 0.25,mu);
exercise_boundary(:,5) = AmericanPut(N, S0, K, T, r, 0.3,mu);

figure(2);
plot(t, exercise_boundary);
ylim([6 10]);
xlabel('time t')
ylabel('S_t')
leg = legend('\sigma=0.1','\sigma=0.15','\sigma=0.2','\sigma=0.25','\sigma=0.3');
title("Exercise Boundary for r=0.02 with various \sigma");

exercise_boundary=[];
exercise_boundary(:,1) = AmericanPut(N, S0, K, T, 0.01, sigma, mu);
exercise_boundary(:,2) = AmericanPut(N, S0, K, T, 0.015, sigma, mu);
exercise_boundary(:,3) = AmericanPut(N, S0, K, T, 0.02, sigma, mu);
exercise_boundary(:,4) = AmericanPut(N, S0, K, T, 0.025, sigma, mu);
exercise_boundary(:,5) = AmericanPut(N, S0, K, T, 0.03, sigma, mu);

figure(3);
plot(t, exercise_boundary);
ylim([6 10]);
xlabel('time t');
ylabel('S_t');
leg = legend('r=0.01','r=0.015','r=0.02','r=0.025','r=0.03');
title("Exercise Boundary for \sigma=0.2 with various r");