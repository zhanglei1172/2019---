%% 黄金分割搜索
clc, clearvars, close all;
dest = 150; %期望100MPa
global poly data
data = xlsread('../../data/附件3-弹性模量与压力.xlsx');
poly = polyfit(data(:, 1), data(:, 2), 4);
P_start = find(data(:, 1)==100);
len = length(data) - P_start+1;
data(P_start, 2);
data_rho = zeros(length(data), 1);
data_rho(P_start) = 0.850;
for i = 1:len-1
    data_rho(P_start+i) = data_rho(P_start+i-1) + ...
        (data_rho(P_start+i-1)*0.5/data(P_start+i-1, 2));
end

for i = -1:-1:1-P_start
    data_rho(P_start+i) = data_rho(P_start+i+1) + ...
        (data_rho(P_start+i+1)*(-0.5)/data(P_start+i+1, 2));
end
data = [data, data_rho];
RHO = data(data(:, 1) == 160, 3); % 高压油泵处密度
%% 设置参数
global gaoYaPipe_v
gaoYaPipe_length = 500; % 500 mm长
gaoYaPipe_diameter = 10;
gaoYaPipe_v = gaoYaPipe_length * (pi*(gaoYaPipe_diameter/2)^2); % 体积
A_diameter = 1.4;
% 单向阀

% 喷油器

% 公式参数
C = 0.85;
A = pi*(A_diameter/2)^2;
Q_func = @(delta_P)C*A*sqrt(2*delta_P/RHO);

delta_t = 0.1;
period = 3e3; % 1s  超参%%%
%% 参数设置
global P0 rho weight
weight = logspace(1, 10, period/delta_t+1); %超参%%%
[P0, rho]  = deal(zeros(2+period/delta_t, 1));
T_change_eps = 5e-4;
T_serach_min = 0.1;
T_serach_max = 1.5;
search_range = T_serach_max - T_serach_min;
%%
T_c = T_serach_min + (1-0.618) * search_range;
T_d = T_serach_min + 0.618 * search_range;
F_c = funEval_tmp(T_c, delta_t, period, Q_func, RHO, dest);
F_d = funEval_tmp(T_d, delta_t, period, Q_func, RHO, dest);
tic
record = [];
res_record = [];
while search_range > T_change_eps
   if(F_c < F_d)
       record = [record T_serach_max T_serach_min];
       res_record = [res_record F_c F_d ];
        T_serach_max = T_d;
        T_d = T_c;
        F_d = F_c;
        search_range = T_serach_max-T_serach_min;
        T_c = T_serach_min+0.382*search_range;
        F_c = funEval_tmp(T_c, delta_t, period, Q_func, RHO, dest);
    else 
%         T_serach_max = T_serach_max;
        T_serach_min = T_c;
        T_c = T_d;
        F_c = F_d;
        search_range = T_serach_max-T_serach_min;
        T_d = T_serach_min+0.618*search_range;
        F_d = funEval_tmp(T_d, delta_t, period, Q_func, RHO, dest);
    end 
end
T_ans = (T_serach_max+T_serach_min)/2
% res = funEval_tmp(T_ans, delta_t, period, Q_func, RHO, dest);
% figure
% plot(P0(1:end-1))
toc
subplot(121)
plot(record, '-*', 'LineWidth', 1.5)
title('开放时间随迭代次数收敛')
xlabel('迭代次数')
ylabel('T')
subplot(122)
plot(res_record, '-*', 'LineWidth', 1.5)
title('评估函数随迭代次数收敛')
xlabel('迭代次数')
ylabel('fun')