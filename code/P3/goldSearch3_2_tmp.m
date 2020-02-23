%%
clc, clearvars, close all
global data_out data_shape data sai_h beng_rho beng_P P0 base_h beng_h beng_s period_out gaoYaPipe_v
data_shape = xlsread('../../data/附件1-凸轮边缘曲线.xlsx');
data_out = xlsread('../../data/附件2-针阀运动曲线.xlsx');
data_out(end,:) = [];
data = xlsread('../../data/附件3-弹性模量与压力.xlsx');
data_shape(end+1, :) = [3.14*2 data_shape(1, 2)]; % 补上2pi 便于计算
%%
omiga = 0.055; % rad/ms
time_delay = 50; % 喷嘴2的相对喷嘴1延后开启工作时间[0, 100]
threshold = 105; % 减压阀阈值

time_delta = 100 - time_delay;
%% 计算rho和P之间的关系
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
%% 设置参数
gaoYaPipe_length = 500; % 500 mm长
gaoYaPipe_diameter = 10;
gaoYaPipe_v = gaoYaPipe_length * (pi*(gaoYaPipe_diameter/2)^2); % 体积
A_diameter = 1.4;
% P_A = 160; % 160 MPa
% 单向阀


% 喷油器

% 油泵
beng_diameter = 5;
beng_s = pi*(beng_diameter/2)^2; % 底面积
beng_up_h = max(data_shape(:, 2)) - min(data_shape(:, 2));
beng_h = 20/beng_s + beng_up_h;
base_h = min(data_shape(:, 2)); % 下止点对应的凸轮径
base_angle = data_shape(data_shape(:, 2)==base_h, 1); % 下止点对应的凸轮角度
% beng_P(1) = 
% 公式参数
C = 0.85;
A = pi*(A_diameter/2)^2;
Q_func = @(delta_P, RHO)C*A*sqrt(2*delta_P/RHO); % RHO高压处的密度
global Q_function
Q_function = @(delta_P, RHO, S)C*S*sqrt(2*delta_P/RHO);
%%
delta_t = 0.1;  % 0.1ms
period_out = 100; % 100ms
% period = period_in * period_out;
period = 3e3; % 1s
% 初始化
global rho
[P0, rho, sai_h, beng_rho] = deal(zeros(2+period/delta_t, 1));
[beng_P]  = deal(zeros(2+period/delta_t, 1));
P0(1) = 100; % 100MPa
rho(1) = 0.850;
beng_P(1) = 0.5;
beng_rho(1) = data(data(:, 1) == beng_P(1), 3); % 高压油泵处初始密度
sai_h(1) = 0;  % 初始塞子偏移
tic

% weight = logspace(1, 10, period/delta_t+1); %超参%%%
omiga_change_eps = 1e-4;
omiga_serach_min = 0.04;
omiga_serach_max = 0.15;
search_range_omiga = omiga_serach_max - omiga_serach_min;

time_delay_change_eps = 1e-1;
time_delay_serach_min = 50;
time_delay_serach_max = 50;
search_range_time_delay = time_delay_serach_max - time_delay_serach_min;

threshold_change_eps = 1e-2;
threshold_serach_min = 101;
threshold_serach_max = 105;
search_range_threshold = threshold_serach_max - threshold_serach_min;

param_serach_min = [omiga_serach_min, time_delay_serach_min, threshold_serach_min];
param_serach_max = [omiga_serach_max, time_delay_serach_max, threshold_serach_max];
param_change_eps = [omiga_change_eps, time_delay_change_eps, threshold_change_eps];
search_range_param = [search_range_omiga, search_range_time_delay, search_range_threshold];
%%
cycle = 2;
T_c = param_serach_min + (1-0.618) * search_range_param;
T_d = param_serach_min + 0.618 * search_range_param;
F_c = funEval2(T_c, delta_t, period, Q_func);
F_d = funEval2(T_d, delta_t, period, Q_func);
tic
record = [];
res_record = [];
while any(search_range_param > param_change_eps)
    cycle = mod(cycle, 3)+1;
    if search_range_param(cycle) < param_change_eps(cycle)
        continue
    end
   record = [record;[param_serach_max param_serach_min]];
   res_record = [res_record F_c F_d ];
   if(F_c < F_d)
        param_serach_max(cycle) = T_d(cycle);
        T_d(cycle) = T_c(cycle);
        F_d = F_c;
        search_range_param(cycle) = param_serach_max(cycle)-param_serach_min(cycle);
        T_c(cycle) = param_serach_min(cycle)+0.382*search_range_param(cycle);
        F_c = funEval2(T_c, delta_t, period, Q_func);
    else 
%         T_serach_max = T_serach_max;
        param_serach_min(cycle) = T_c(cycle);
        T_c(cycle) = T_d(cycle);
        F_c = F_d;
        search_range_param(cycle) = param_serach_max(cycle)-param_serach_min(cycle);
        T_d(cycle) = param_serach_min(cycle)+0.618*search_range_param(cycle);
        F_d = funEval2(T_d, delta_t, period, Q_func);
   end
   if abs(F_c - F_d) < 1e-2
       search_range_param(cycle) = 0;
   end
end
T_ans_add = (param_serach_max+param_serach_min)/2
% res = funEval(T_ans_add, delta_t, period, Q_func, RHO, dest);
toc
subplot(121)
plot(record, '-*', 'LineWidth', 1.5)
title('\omega随迭代次数收敛')
xlabel('迭代次数')
ylabel('\omega')
subplot(122)
plot(res_record, '-*', 'LineWidth', 1.5)
title('评估函数随迭代次数收敛')
xlabel('迭代次数')
ylabel('fun')