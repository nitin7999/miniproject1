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

    y = X * F * h + n;
    H = X*F;

    h_est = inv(H' * H) * H' * y;
    
    error_sum = error_sum + (h_est - h)' * (h_est - h);
    h_sum = h_sum + h_est;
    MSE_th_sum = MSE_th_sum + trace( sigma^2 * inv(H' * H));
end

average_error = error_sum/10000;
h_avg = h_sum/10000;
MSE_th = MSE_th_sum/10000;

plot(abs(h - h_est));


