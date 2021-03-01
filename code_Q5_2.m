clc
clear all;
close all;
N = 512; % number of symbols
Nb_s = 2; % number of bits per symbol
Nb = N * Nb_s; % number of bits
L = 32;

lambda = 0.2;
sigma = .1;

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

h = p .* ((2^-0.5) * (randn(L,1) + 1j * randn(L,1))) ./(norm(p));

h_sparsity_index = randperm(32,6);
h_sparsity_index = sort(h_sparsity_index);
h_sparsity = zeros(L,1);
h_sparsity(h_sparsity_index) = h(h_sparsity_index);

y = X * F * h_sparsity + n;
H = X*F;

r = y;
h_sparsity_index_est = [];

for index = 1:6
    
    [temp temp_I] = max(H'*r);
    h_sparsity_index_est = [h_sparsity_index_est temp_I];
    for index1 = 1:length(h_sparsity_index_est)
        A_temp(:,index1) =    H(:,h_sparsity_index_est(index1));
    end
    P = A_temp * inv(A_temp' * A_temp) * A_temp';
    
    r = (eye(N) -  P) * y;
    
end
    
b = zeros(L - 6 , 1);
sample_space = 1:32;
A_indices = setdiff(sample_space,h_sparsity_index_est);
A = zeros(L - 6 , L);
for i = 1:26
    A(i,A_indices(i)) = 1;
end


h_est_temp = inv(H' * H) * H' * y;
h_est_sparcity = inv(H' * H) * H' * y - inv(H' * H) * A' * inv(A * inv(H' * H) * A') * A * inv(H' * H) * H' * y;


stem(abs(h_sparsity)); figure; stem(abs(h_est_sparcity));figure; stem(abs(h_est_sparcity - h_sparsity));
