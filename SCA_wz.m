function [W,Z,SEE] = SCA_wz(K,N_T,P_max,P_c,noise,xi,bar_g,bar_g_eve,bar_a,bar_b,bar_f,bar_beta)
%UNTITLED3 此处提供此函数的摘要
%   此处提供详细说明
num = 1;
opt(1) = inf;
while 1
    lambda = 0;
    while 1
    
        cvx_begin quiet%防止模型在求解时生成任何屏幕输出
        cvx_solver sdpt3%求解器
        variable W(N_T,N_T,K) complex semidefinite%变量W 复杂半正定矩阵
        variable Z(N_T,N_T) complex semidefinite%变量Z 复杂半正定矩阵
        variable alpha_v(K)%users效率的仿射received information rate
        variable beta_v(K)%eavesdropper效率的仿射information leakage rate
        variable a(K)%2的alpha次方-1的e次开方
        variable b(K)%
        variable c(K)
        variable d(K)
        variable f(K)%有待商榷
        expression total_power%表达式 总能量
        expression sum_secure_rate%表达式 总保密能效
        expression sum_matrix%表达式 总矩阵
        sum_secure_rate = 0;%初始化
        sum_matrix = zeros(N_T,N_T);%初始化
        for i = 1:K
            sum_matrix = sum_matrix + W(:,:,i);%1到k个用户的波束赋形向量总和
            sum_secure_rate = sum_secure_rate + (alpha_v(i) - beta_v(i));%1到k个用户的保密能效总和
        end
        total_power = trace(sum_matrix + Z) + P_c;%总能量=（1到k个用户的波束赋形向量总和+噪声的波束赋形向量）的迹+电路模块消耗的常量（constant power consumption）
        maximize sum_secure_rate - lambda * total_power%由lambda=sum_secure_rate/total_power变换得到
        subject to %约束条件
        for i = 1:K%PA3
            alpha_v(i) >= xi;%xi是QoS requirement用户所需信息速率
            bar_g(i,:)*W(:,:,i)*bar_g(i,:)' >= exp(a(i) + b(i));%
            bar_g_eve*W(:,:,i)*bar_g_eve' <= exp(bar_f(i)) * (1+ f(i) - bar_f(i));%
            c(i) + d(i) >= f(i);%
            exp(bar_a(i)) * (1 + a(i) -bar_a(i)) >= 2^alpha_v(i) - 1;%
            exp(bar_b(i)) * (1 + b(i) -bar_b(i)) >= bar_g(i,:)*(sum_matrix + Z - W(:,:,i))*bar_g(i,:)' + noise;%
            (2^bar_beta(i)) + (2^bar_beta(i))*log(2)*(beta_v(i) - bar_beta(i)) >= exp(c(i)) + 1;%
            bar_g_eve * (sum_matrix + Z - W(:,:,i)) * bar_g_eve' + noise >= exp(d(i));%
        end
        total_power <= P_max;%能量<预算
        cvx_end

        if sum_secure_rate - lambda * total_power < 10^-4%
            break;
        end

        lambda = sum_secure_rate/total_power;
    
    end

    for i = 1:K
        r_temp(i) = real(log2(1 + bar_g(i,:)*W(:,:,i)*bar_g(i,:)'/(bar_g(i,:)*(sum_matrix + Z - W(:,:,i))*bar_g(i,:)' + noise)));%
        r_eve_temp(i) = real(log2(1 + bar_g_eve*W(:,:,i)*bar_g_eve'/(bar_g_eve*(sum_matrix + Z - W(:,:,i))*bar_g_eve' + noise)));%
    end
        
    for i = 1:K%（31）
        bar_a(i) = log(2^r_temp(i) - 1);%
        bar_b(i) = real(log(bar_g(i,:)*(sum_matrix-W(:,:,i)+Z)*bar_g(i,:)' + noise));%
        bar_c(i) = log(2^r_eve_temp(i) - 1);%
        bar_d(i) = real(log(bar_g_eve*(sum_matrix + Z - W(:,:,i))*bar_g_eve' + noise));%
        bar_f(i) = bar_c(i) + bar_d(i);%
        bar_beta(i) = r_eve_temp(i);%
    end

    num = num + 1;%
    SEE_temp(num) = sum(r_temp - r_eve_temp)/total_power%第几个SEE
    
    if (SEE_temp(num) - SEE_temp(num-1)) < 10^-4%到达最大值后停下
        break;
    end

end

SEE_seq = SEE_temp(2:length(SEE_temp));
SEE = SEE_temp(length(SEE_temp));
end