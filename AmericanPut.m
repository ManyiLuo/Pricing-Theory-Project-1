function [exercise_boundary,optionPrice] = AmericanPut(N,S0,K,T,r,sigma,mu)
% tree parameters
dt = T/N;         % size of timesteps
t = [0:dt:T];
u = exp(r*dt+sigma*sqrt(dt));  % stock price up ratio
d = exp(r*dt-sigma*sqrt(dt));                  % stock price down ratio
%dscnt = exp(-r*dt);       % discount rate
qu = (1-exp(-sigma*sqrt(dt)))/(exp(sigma*sqrt(dt))-exp(-sigma*sqrt(dt))); % probability (up)
qd = 1 - qu;   % probability (down)

% First we need an empty matrix big enough to hold all the calculated prices
priceTree = nan(N + 1, N + 1); % matrix filled with NaN == missing values
% Note that this tree is approx twice as big as needed because we only need one side of diagonal
% We use the top diagonal for efficiency

% Initialize with the current stock price, then loop through all the steps
priceTree(1, 1) = S0;
for ii = 2:(N+1)  % add the endfor now
  % vector calculation of all the up steps (show how the indexing works on line)
  priceTree(1:(ii-1), ii) = priceTree(1:(ii-1), ii-1) * u;

  % The diagonal will hold the series of realizations that is always down
  % this is a scalar calculation
  priceTree(ii, ii) = priceTree((ii-1), (ii-1)) * d;
end

% Now we can calculate the value of the option at each node
% We need a matrix to hold the option values that is the same size as the price tree
% Note that size returns a vector of dimensions [r, c] which is why it can be passed as dimensions to nan
optionValueTree = nan(size(priceTree));

% First we calculate the end value
optionValueTree(:, end) = max(0, K - priceTree(:, end));

%exercise boundary
%exercise_boundary = zeros(1,N+1);
%exercise_boundary(N+1) = priceTree(find(optionValueTree(:, end) == K - priceTree(:, end), 1, 'first'), N+1);

discountRate = exp(-r * dt);

exercise_boundary = nan(1,N+1);
exercise_boundary(N+1) = priceTree(find(optionValueTree(:, end) == K - priceTree(:, end), 1, 'first'), N+1);

for ii = N:-1:1
    % If do not exercise, then the value is
    optionValueTree(1:ii, ii) = discountRate * (qu * optionValueTree(1:ii, ii+1) + qd * optionValueTree(2:(ii+1), ii+1));
    % If exercise, then the value is
    exercise_value = K - priceTree(1:ii, ii);
    optionValueTree(1:ii, ii) = max(optionValueTree(1:ii,ii), exercise_value);
    
    jj = find(optionValueTree(1:ii, ii) == exercise_value, 1, 'first');
    if isempty(jj) == 0
        exercise_boundary(ii) = priceTree(jj,ii);
    end
    
end

optionPrice = optionValueTree(1, 1)

end