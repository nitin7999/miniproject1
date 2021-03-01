clc
clear all;
close all;
N = 512; % number of symbols
Nb_s = 2; % number of bits per symbol
Nb = N * Nb_s; % number of bits
L = 32;

lambda = 0.2;
sigma = 1;

j = 1:L;

p = exp( -lambda * (j - 1).');

h = p .* ((2^-0.5) * (randn(L,1) + 1j * randn(L,1))) ./(norm(p));

h_sparsity_index = randperm(32,6);
h_sparsity_index = sort(h_sparsity_index);
h_sparsity = zeros(L,1);
h_sparsity(h_sparsity_index) = h(h_sparsity_index);

b = zeros(L - 6 , 1);
sample_space = 1:32;
A_indices = setdiff(sample_space,h_sparsity_index);
A = zeros(L - 6 , L);
for i = 1:26
    A(i,A_indices(i)) = 1;
end

error_sum = 0;
h_sum = 0;
MSE_th_sum = 0;
i = 1:N; 

for index = 1:10000

    bits = round(rand(Nb,1)); %random bits generation.
    H = eye(N);

    c_index = 2*bits(1:2:end) + bits(2:2:end) + 1;

    C = 2^(-0.5) * [1+1j , 1-1j , -1+1j, -1-1j];

    S = C(c_index); % symbol set
    X = S.*H;

    
    temp1 = repmat( exp(2*pi * 1i * (i.' - 1) / 512),1,L);
    temp2 = repmat(j,N,1);

    F = temp1.^(temp2-1);

    n = sigma * (2^-0.5) * (randn(N,1) + 1j * randn(N,1));

    y = X * F * h_sparsity + n;
    H = X*F;

    h_est_sparcity = inv(H' * H) * H' * y - inv(H' * H) * A' * inv(A * inv(H' * H) * A') * A * inv(H' * H) * H' * y;

    error_sum = error_sum + (h_est_sparcity - h_sparsity)' * (h_est_sparcity - h_sparsity);
    h_sum = h_sum + h_est_sparcity;
    MSE_th_sum = MSE_th_sum + trace( sigma^2 * (inv(H' * H) - inv(H' * H) * A' * inv( A * inv(H' * H) * A') * A * inv(H' * H)) );

end

average_error = error_sum/10000;
h_avg = h_sum/10000;
MSE_th = MSE_th_sum/10000;
