% Initialize Parameters
N = 10; % Number of macro channels
K = 5; % Number of small base stations (SBS)
M = 15; % Number of secondary users (SUs)
P_mbs = 20; % Transmission power of MBS in Watt
P_sbs = 1; % Maximum transmission power of SBS in Watt
sigma = 0.0001; % Noise power
lambda = 0.5; % Weight parameter
alpha = 0.01; % Price incrementing factor
prob_mbs_on = 0.75; % Probability of MBS being ON

% Initialize random channel gains
h_km = rand(K, M); % Channel gain between SBS k and SU m
h_nk = rand(N, K); % Channel gain between MBS using channel n and SBS k's user
h_ni = rand(N, M); % Channel gain between MBS n and its MPU i
h_ki = rand(K, M); % Interference channel gain between SBS k and MPU i

% Initialize utility matrices
U_joint = zeros(K, N);
C_km = zeros(K, M);

% Initialize price vector
p_n = ones(1, N) * 0.1; % Initial price of channels

% SINR calculation
SINR = @(x_km, p_k, h_km) (x_km * p_k * abs(h_km)^2) / sigma;

% Channel capacity calculation using Shannon's formula
channel_capacity = @(SINR) log2(1 + SINR);

% Utility calculations
for k = 1:K
    for n = 1:N
        % Calculate channel capacities for SBS
        for m = 1:M
            x_km = 1; % Assume SU is associated with SBS k
            SINR_val = SINR(x_km, P_sbs, h_km(k, m));
            C_km(k, m) = channel_capacity(SINR_val);
        end
        
        % Utility with respect to SBS channels
        R_k_n_SBS = prob_mbs_on * channel_capacity((P_mbs * abs(h_nk(n, k))^2) / (sigma + P_mbs * abs(h_ni(n, k))^2));
        U_k_n_SBS = (1 / M) * sum(R_k_n_SBS);
        
        % Utility with respect to MBS channels
        R_k_n_MBS = prob_mbs_on * channel_capacity((P_mbs * abs(h_ni(n, k))^2) / (sigma + P_sbs * abs(h_ki(k, n))^2));
        U_k_n_MBS = R_k_n_MBS;
        
        % Joint utility
        U_joint(k, n) = lambda * U_k_n_SBS + (1 - lambda) * U_k_n_MBS;
        
        fprintf('The Joint Utility of the primary and the secondary system is:',U_joint);
    end
end

% Auction-based algorithm
D_k = cell(1, K); % Demand sets
D_ex = []; % Excess demand set

count = 0;
while true
    count = count + 1;
    
    % Stage I: Demand set calculation
    for k = 1:K
        monetary_gain = U_joint(k, :) - p_n;
        [~, sorted_indices] = sort(monetary_gain, 'descend');
        D_k{k} = sorted_indices(1:min(M, numel(sorted_indices))); % Ensure not to exceed the number of elements
    end
    
    % Stage II: Excess demand set calculation
    D_ex = [];
    for n = 1:N
        demand_count = 0;
        for k = 1:K
            if any(D_k{k} == n)
                demand_count = demand_count + 1;
            end
        end
        if demand_count > 1
            D_ex = [D_ex, n];
        end
    end
    
    % Break if there is no excess demand
    if isempty(D_ex)
        break;
    end
    
    % Stage III: Walrasian Equilibrium implementation
    for n = D_ex
        p_n(n) = p_n(n) + alpha;
    end
end

% Calculate sum rates for SBS and MBS
sum_rate_SBS = zeros(1, K);
sum_rate_MBS = zeros(1, N);

for k = 1:K
    for idx = 1:numel(D_k{k})
        n = D_k{k}(idx);
        for m = 1:M
            x_km = 1; % Assume SU is associated with SBS k
            SINR_val = SINR(x_km, P_sbs, h_km(k, m));
            sum_rate_SBS(k) = sum_rate_SBS(k) + channel_capacity(SINR_val);
        end
    end
end

for n = 1:N
    for k = 1:K
        if any(D_k{k} == n)
            sum_rate_MBS(n) = sum_rate_MBS(n) + prob_mbs_on * channel_capacity((P_mbs * abs(h_nk(n, k))^2) / (sigma + P_mbs * abs(h_ni(n, k))^2));
        end
    end
end

% Plotting
figure;
scatter(sum(sum_rate_SBS), sum(sum_rate_MBS), 'filled');
xlabel('SBS Sum Rate');
ylabel('MBS Sum Rate');
title('Sum Rate of SBS vs MBS');
legend('Walrasian Equilibrium');
grid on;
