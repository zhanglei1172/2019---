%% 找到的最优omiga下找最优的相对后延时间
clc, clearvars, close all
global  data_shape data sai_h beng_rho P0 base_h beng_h beng_s period_out gaoYaPipe_v
data_shape = xlsread('../../data/附件1-凸轮边缘曲线.xlsx');
data_out = xlsread('../../data/附件2-针阀运动曲线.xlsx');
data_out(end,:) = [];
data = xlsread('../../data/附件3-弹性模量与压力.xlsx');
data_shape(end+1, :) = [3.14*2 data_shape(1, 2)]; % 补上2pi 便于计算
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
profile on
delta_t = 0.1;  % 0.1ms
period_out = 100; % 100ms
% period = period_in * period_out;
period = 3e3; % 1s
% 初始化
global rho delta_theta
[P0, rho, sai_h, beng_rho] = deal(zeros(2+period/delta_t, 1));
% [beng_P]  = deal(zeros(2+period/delta_t, 1));
P0(1) = 100; % 100MPa
rho(1) = 0.850;
beng_P(1) = 0.5;
beng_rho(1) = data(data(:, 1) == beng_P(1), 3); % 高压油泵处初始密度
sai_h(1) = 0;  % 初始塞子偏移
omiga = 0.05545;
delta_theta = delta_t * omiga; % 每次增加的角度
tic
T_change_eps = 1e-3;
%%
T_serach_min = 30;
T_serach_max = 70;
search_range = T_serach_max - T_serach_min;
T_c = T_serach_min + (1-0.618) * search_range;
T_d = T_serach_min + 0.618 * search_range;
F_c = funEvalT(delta_t, period, Q_func, T_c,  data_out);
F_d = funEvalT(delta_t, period, Q_func, T_d, data_out);
tic
record = [];
res_record = [];
while search_range > T_change_eps
   record = [record T_serach_max T_serach_min];
   res_record = [res_record F_c F_d ];
   if(F_c < F_d)
        T_serach_max = T_d;
        T_d = T_c;
        F_d = F_c;
        search_range = T_serach_max-T_serach_min;
        T_c = T_serach_min+0.382*search_range;
        F_c = funEvalT(delta_t, period, Q_func, T_c, data_out);
    else 
%         T_serach_max = T_serach_max;
        T_serach_min = T_c;
        T_c = T_d;
        F_c = F_d;
        search_range = T_serach_max-T_serach_min;
        T_d = T_serach_min+0.618*search_range;
        F_d = funEvalT(delta_t, period, Q_func, T_d,  data_out);
   end
   
end

% res = funEval(T_ans_add, delta_t, period, Q_func, RHO, dest);

% profile viewer
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
disp('最佳相对后沿时间：')
T_ans_add = (T_serach_max+T_serach_min)/2