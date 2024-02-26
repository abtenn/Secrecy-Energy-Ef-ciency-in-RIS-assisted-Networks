function [bar_a,bar_b,bar_f,bar_beta] = SCA_initialization(K,N_T,P_max,P_c,noise,xi,bar_g,bar_g_eve)
%UNTITLED 此处提供此函数的摘要
%   此处提供详细说明
cvx_begin quiet
cvx_solver sdpt3
variable W(N_T,N_T,K) complex semidefinite
variable Z(N_T,N_T) complex semidefinite
expression total_power_temp
expression sum_matrix_temp
sum_matrix_temp = zeros(N_T,N_T);
for i = 1:K
    sum_matrix_temp = sum_matrix_temp + W(:,:,i);
end
total_power_temp = trace(Z);
for i = 1:K
    total_power_temp = total_power_temp + trace(W(:,:,i));
end
subject to %约束条件
    for i = 1:K
        bar_g(i,:)*W(:,:,i)*bar_g(i,:)' >= (2^xi - 1) * (bar_g(i,:)*(sum_matrix_temp + Z - W(:,:,i))*bar_g(i,:)' + noise);
    end
    total_power_temp + P_c <= P_max;
cvx_end

for i = 1:K
    r_temp(i) = real(log2(1 + bar_g(i,:)*W(:,:,i)*bar_g(i,:)'/(bar_g(i,:)*(sum_matrix_temp + Z - W(:,:,i))*bar_g(i,:)' + noise)));
    r_eve_temp(i) = real(log2(1 + bar_g_eve*W(:,:,i)*bar_g_eve'/(bar_g_eve*(sum_matrix_temp + Z - W(:,:,i))*bar_g_eve' + noise)));  
end

for i = 1:K
    bar_a(i) = log(2^r_temp(i) - 1);
    bar_b(i) = real(log(bar_g(i,:)*(sum_matrix_temp-W(:,:,i)+Z)*bar_g(i,:)' + noise));
    bar_c(i) = log(2^r_eve_temp(i) - 1);
    bar_d(i) = real(log(bar_g_eve*(sum_matrix_temp + Z - W(:,:,i))*bar_g_eve' + noise));
    bar_f(i) = bar_c(i) + bar_d(i);
    bar_beta(i) = r_eve_temp(i);
end
end