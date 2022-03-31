% computation parameter
N = 5000;  % number of timesteps

% option parameters
S0 = 10;      % initial stock price
K = 10;       % strike price
T = 1;       % expiry time
r = 0.02;      % interest rate
sigma = 0.1;  % volatility
mu = 0.05; %drift

% tree parameters
dt = T/N;         % size of timesteps
t = [0:dt:T];
u = exp(r*dt+sigma*sqrt(dt));  % stock price up ratio
d = exp(r*dt-sigma*sqrt(dt));                  % stock price down ratio
%dscnt = exp(-r*dt);       % discount rate
pu = 0.5 * (1 + (mu - r - 0.5 * sigma^2)/sigma * sqrt(dt)); % probability (up)
pd = 1 - pu;   % probability (down)

sim = 10000;
exercise_time = NaN(1, sim);
exercise_profit = NaN(1, sim);
exercise_boundary = AmericanPut(N,S0,K,T,r,sigma,mu);
for i = 1:sim
    S = NaN(1, N+1);
    S(1) = S0;

    for j = 1:N
        indicator = rand;
            if indicator <= pu
                y = 1;
            else
                y = -1;
            end
        
        S(j+1) = S(j) * exp(r * dt + sigma * sqrt(dt) * y);

        if S(j+1) <= exercise_boundary(j)
            exercise_time(i) = j * dt;
            exercise_profit(i) = max(0, K - S(j)) * exp(-r * j * dt);
            break;
        end

        if j == N
            exercise_time(i) = 1;
            exercise_profit(i) = max(0, K - S(N+1)) * exp(- r * T);
        end
    end
end

exercise_profit = exercise_profit - 0.3225;
exercise_time = exercise_time(exercise_time ~=1);
exercise_profit = exercise_profit(exercise_profit ~=-0.3225);

%p1 = ksdensity(exercise_profit);
%p2 = ksdensity(exercise_time);

[g,yi] = ksdensity(exercise_time);
mean1 = mean([g,yi]);
sd1 = std([g,yi]);
q1 = quantile([g,yi],[0.25 0.50 0.75])
plot(yi,g,'Linewidth',1);
hold on

% computation parameter
N = 5000;  % number of timesteps

% option parameters
S0 = 10;      % initial stock price
K = 10;       % strike price
T = 1;       % expiry time
r = 0.02;      % interest rate
sigma = 0.15;  % volatility
mu = 0.05; %drift

% tree parameters
dt = T/N;         % size of timesteps
t = [0:dt:T];
u = exp(r*dt+sigma*sqrt(dt));  % stock price up ratio
d = exp(r*dt-sigma*sqrt(dt));                  % stock price down ratio
%dscnt = exp(-r*dt);       % discount rate
pu = 0.5 * (1 + (mu - r - 0.5 * sigma^2)/sigma * sqrt(dt)); % probability (up)
pd = 1 - pu;   % probability (down)

sim = 10000;
exercise_time = NaN(1, sim);
exercise_profit = NaN(1, sim);
exercise_boundary = AmericanPut(N,S0,K,T,r,sigma,mu);
for i = 1:sim
    S = NaN(1, N+1);
    S(1) = S0;

    for j = 1:N
        indicator = rand;
            if indicator <= pu
                y = 1;
            else
                y = -1;
            end
        
        S(j+1) = S(j) * exp(r * dt + sigma * sqrt(dt) * y);

        if S(j+1) <= exercise_boundary(j)
            exercise_time(i) = j * dt;
            exercise_profit(i) = max(0, K - S(j)) * exp(-r * j * dt);
            break;
        end

        if j == N
            exercise_time(i) = 1;
            exercise_profit(i) = max(0, K - S(N+1)) * exp(- r * T);
        end
    end
end

exercise_profit = exercise_profit - 0.5162;
exercise_time = exercise_time(exercise_time ~=1);
exercise_profit = exercise_profit(exercise_profit ~=-0.5162);

[g,yi] = ksdensity(exercise_time);
mean2 = mean([g,yi]);
sd2 = std([g,yi]);
q2 = quantile([g,yi],[0.25 0.50 0.75])
plot(yi,g,'Linewidth',1);
hold on

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
pu = 0.5 * (1 + (mu - r - 0.5 * sigma^2)/sigma * sqrt(dt)); % probability (up)
pd = 1 - pu;   % probability (down)

sim = 10000;
exercise_time = NaN(1, sim);
exercise_profit = NaN(1, sim);
exercise_boundary = AmericanPut(N,S0,K,T,r,sigma,mu);
for i = 1:sim
    S = NaN(1, N+1);
    S(1) = S0;

    for j = 1:N
        indicator = rand;
            if indicator <= pu
                y = 1;
            else
                y = -1;
            end
        
        S(j+1) = S(j) * exp(r * dt + sigma * sqrt(dt) * y);

        if S(j+1) <= exercise_boundary(j)
            exercise_time(i) = j * dt;
            exercise_profit(i) = max(0, K - S(j)) * exp(-r * j * dt);
            break;
        end

        if j == N
            exercise_time(i) = 1;
            exercise_profit(i) = max(0, K - S(N+1)) * exp(- r * T);
        end
    end
end

exercise_profit = exercise_profit - 0.7111;
exercise_time = exercise_time(exercise_time ~=1);
exercise_profit = exercise_profit(exercise_profit ~=-0.7111);

[g,yi] = ksdensity(exercise_time);
mean3 = mean([g,yi]);
sd3 = std([g,yi]);
q3 = quantile([g,yi],[0.25 0.50 0.75])
plot(yi,g,'Linewidth',1);
hold on


% computation parameter
N = 5000;  % number of timesteps

% option parameters
S0 = 10;      % initial stock price
K = 10;       % strike price
T = 1;       % expiry time
r = 0.02;      % interest rate
sigma = 0.25;  % volatility
mu = 0.05; %drift

% tree parameters
dt = T/N;         % size of timesteps
t = [0:dt:T];
u = exp(r*dt+sigma*sqrt(dt));  % stock price up ratio
d = exp(r*dt-sigma*sqrt(dt));                  % stock price down ratio
%dscnt = exp(-r*dt);       % discount rate
pu = 0.5 * (1 + (mu - r - 0.5 * sigma^2)/sigma * sqrt(dt)); % probability (up)
pd = 1 - pu;   % probability (down)

sim = 10000;
exercise_time = NaN(1, sim);
exercise_profit = NaN(1, sim);
exercise_boundary = AmericanPut(N,S0,K,T,r,sigma,mu);
for i = 1:sim
    S = NaN(1, N+1);
    S(1) = S0;

    for j = 1:N
        indicator = rand;
            if indicator <= pu
                y = 1;
            else
                y = -1;
            end
        
        S(j+1) = S(j) * exp(r * dt + sigma * sqrt(dt) * y);

        if S(j+1) <= exercise_boundary(j)
            exercise_time(i) = j * dt;
            exercise_profit(i) = max(0, K - S(j)) * exp(-r * j * dt);
            break;
        end

        if j == N
            exercise_time(i) = 1;
            exercise_profit(i) = max(0, K - S(N+1)) * exp(- r * T);
        end
    end
end

exercise_profit = exercise_profit - 0.9063;
exercise_time = exercise_time(exercise_time ~=1);
exercise_profit = exercise_profit(exercise_profit ~=-0.9063);

[g,yi] = ksdensity(exercise_time);
mean4 = mean([g,yi]);
sd4 = std([g,yi]);
q4 = quantile([g,yi],[0.25 0.50 0.75]);
plot(yi,g,'Linewidth',1);
hold on


% computation parameter
N = 5000;  % number of timesteps

% option parameters
S0 = 10;      % initial stock price
K = 10;       % strike price
T = 1;       % expiry time
r = 0.02;      % interest rate
sigma = 0.3;  % volatility
mu = 0.05; %drift

% tree parameters
dt = T/N;         % size of timesteps
t = [0:dt:T];
u = exp(r*dt+sigma*sqrt(dt));  % stock price up ratio
d = exp(r*dt-sigma*sqrt(dt));                  % stock price down ratio
%dscnt = exp(-r*dt);       % discount rate
pu = 0.5 * (1 + (mu - r - 0.5 * sigma^2)/sigma * sqrt(dt)); % probability (up)
pd = 1 - pu;   % probability (down)

sim = 10000;
exercise_time = NaN(1, sim);
exercise_profit = NaN(1, sim);
exercise_boundary = AmericanPut(N,S0,K,T,r,sigma,mu);
for i = 1:sim
    S = NaN(1, N+1);
    S(1) = S0;

    for j = 1:N
        indicator = rand;
            if indicator <= pu
                y = 1;
            else
                y = -1;
            end
        
        S(j+1) = S(j) * exp(r * dt + sigma * sqrt(dt) * y);

        if S(j+1) <= exercise_boundary(j)
            exercise_time(i) = j * dt;
            exercise_profit(i) = max(0, K - S(j)) * exp(-r * j * dt);
            break;
        end

        if j == N
            exercise_time(i) = 1;
            exercise_profit(i) = max(0, K - S(N+1)) * exp(- r * T);
        end
    end
end

exercise_profit = exercise_profit - 1.1014;
exercise_time = exercise_time(exercise_time ~=1);
exercise_profit = exercise_profit(exercise_profit ~=-1.1014);

[g,yi] = ksdensity(exercise_time);
mean5 = mean([g,yi]);
sd5 = std([g,yi]);
q5 = quantile([g,yi],[0.25 0.50 0.75])
plot(yi,g,'Linewidth',1);
hold off
xlabel('Exercise Time')
ylabel('Probability Density')
leg = legend('\sigma=0.1','\sigma=0.15','\sigma=0.2','\sigma=0.25','\sigma=0.3');
title("KDE of Exercise Time for r=0.02, \mu = 0.05 and K=10 with various \sigma");

%sigma = {0.1;0.15;0.2;0.25,0.3};
mean = {mean1; mean2; mean3; mean4; mean5};
sd = {sd1; sd2; sd3; sd4; sd5};
q = {q1;q2;q3;q4;q5};
table1 = table(mean,sd,q);