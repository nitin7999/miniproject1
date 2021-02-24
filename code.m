clc
clear all;
close all;
N = 512; % number of symbols
Nb_s = 2; % number of bits per symbol
Nb = N * Nb_s; % number of bits
L = 32;

lambda = 0.2;
sigma = 1;

bits = round(rand(Nb,1)); %random bits generation.
H = eye(N);

c_index = 2*bits(1:2:end) + bits(2:2:end) + 1;

C = 2^(-0.5) * [1+1j , 1-1j , -1+1j, -1-1j];

S = C(c_index); % symbol set
X = S.*H;

i = 1:N; j = 1:L;
temp1 = repmat( exp(2*pi * 1i * (i.' - 1) / 512),1,L);
temp2 = repmat(j,N,1);

F = temp1.^(temp2-1);

n = sigma * (2^-0.5) * (randn(N,1) + 1j * randn(N,1));

p = exp( -lambda * (j - 1).');

h = p .* ((2^-0.5) * (randn(L,1) + 1j * randn(L,1))) ./(abs(p).^2);

y = X * F * h + n;
H = X*F;

h_est = inv(H' * H) * H' * y;



