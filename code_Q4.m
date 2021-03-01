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

b = zeros(3 , 1);

A = zeros(3 , L);
A(1,1:2) = [1 -1];
A(2,3:4) = [1 -1];
A(3,5:6) = [1 -1];

p = exp( -lambda * (j - 1).');

h = p .* ((2^-0.5) * (randn(L,1) + 1j * randn(L,1))) ./(norm(p));
h(2) = h(1);
h(4) = h(3);
h(6) = h(5);

error_sum = 0;
h_sum = 0;
MSE_th_sum = 0;
i = 1:N; 
for index1 = 1:10000

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

    y = X * F * h + n;
    H = X*F;

    h_est = inv(H' * H) * H' * y - inv(H' * H) * A' * inv(A * inv(H' * H) * A') * A * inv(H' * H) * H' * y;

    error_sum = error_sum + (h_est - h)' * (h_est - h);
    h_sum = h_sum + h_est;
    MSE_th_sum = MSE_th_sum + trace( sigma^2 * (inv(H' * H) - inv(H' * H) * A' * inv( A * inv(H' * H) * A') * A * inv(H' * H)) );

end

average_error = error_sum/10000;
h_avg = h_sum/10000;
MSE_th = MSE_th_sum/10000;

stem(abs(h)); figure; stem(abs(h_est));
