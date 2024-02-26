clear;clc;

%% system setting 
K = 4; %number of users
N = 8; %number of RIS
N_T = 8;% number of antennas

P_max = 1;
P_c = 0.1;
noise = 0.01;
xi = 1;%QoS requirement

%% 给H,g,h,g_eve,h_eve和bar_theta赋值

%channel

for i = 1:K
    H = rand(N,N_T) + 1i*rand(N,N_T);%H: Tx-RIS
    g(:,i) = rand(N_T,1) + 1i*rand(N_T,1);%g: Tx-user
    h(:,i) = rand(N,1) + 1i*rand(N,1);%h: RIS-user
end

g_eve =  rand(N_T,1) + 1i*randn(N_T,1);%g: RIS-Eve
h_eve = rand(N,1) + 1i*rand(N,1);%h: Tx-Eve
%theta initialization
for n = 1:N%归一化
    bar_theta(n) =  rand(1,1) + 1i*rand(1,1);
    bar_theta(n) =  bar_theta(n)/sqrt(bar_theta(n)' * bar_theta(n));
end

%% bar_g & bar_g_eve
for i = 1:K
    bar_g(i,:) = (g(:,i)' + h(:,i)'*diag(bar_theta)*H);%全过程
end
bar_g_eve = (g_eve' + h_eve'*diag(bar_theta)*H);

%% initialization of SCA
[bar_a,bar_b,bar_f,bar_beta] = SCA_initialization(K,N_T,P_max,P_c,noise,xi,bar_g,bar_g_eve);

%% optimization of w & z
[W,Z,SEE_w] = SCA_wz(K,N_T,P_max,P_c,noise,xi,bar_g,bar_g_eve,bar_a,bar_b,bar_f,bar_beta);

%% result
SEE_w







