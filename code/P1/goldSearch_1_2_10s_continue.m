%% 黄金分割搜索
clc, clearvars, close all;
dest = 150; %期望100MPa
global poly data second
second = 10e3; % ms
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
period = second; % 1s  超参%%%
period = period + 3e3;
%% 参数设置
global P0 rho weight
% weight = logspace(1, 10, period/delta_t+1); %超参%%%
[P0, rho]  = deal(zeros(2+period/delta_t, 1));
T1_change_eps = 5e-2;
T1_serach_min = 0.28767;
T1_serach_max = 0.28767;
search_range1 = T1_serach_max - T1_serach_min;

T2_change_eps = 1e-8;
T2_serach_min = 1e-8;
T2_serach_max = 1e-5;
search_range2 = T2_serach_max - T2_serach_min;

T_change_eps = [T1_change_eps, T2_change_eps];
T_serach_min = [T1_serach_min, T2_serach_min];
T_serach_max = [T1_serach_max, T2_serach_max];
search_range = T_serach_max - T_serach_min;
%%
cycle = 1;
T_c = T_serach_min + (1-0.618) * search_range;
T_d = T_serach_min + 0.618 * search_range;
F_c = funEval2(T_c, delta_t, period, Q_func, RHO, dest);
F_d = funEval2(T_d, delta_t, period, Q_func, RHO, dest);
tic
record = [];
res_record = [];
while any(search_range > T_change_eps)
    cycle = mod(cycle, 2)+1;
    if search_range(cycle) < T_change_eps(cycle)
        continue
    end
   record = [record;[T_serach_max T_serach_min]];
   res_record = [res_record F_c F_d ];
   if(F_c < F_d)
        T_serach_max(cycle) = T_d(cycle);
        T_d(cycle) = T_c(cycle);
        F_d = F_c;
        search_range(cycle) = T_serach_max(cycle)-T_serach_min(cycle);
        T_c(cycle) = T_serach_min(cycle)+0.382*search_range(cycle);
        F_c = funEval2(T_c, delta_t, period, Q_func, RHO, dest);
    else 
%         T_serach_max = T_serach_max;
        T_serach_min(cycle) = T_c(cycle);
        T_c(cycle) = T_d(cycle);
        F_c = F_d;
        search_range(cycle) = T_serach_max(cycle)-T_serach_min(cycle);
        T_d(cycle) = T_serach_min(cycle)+0.618*search_range(cycle);
        F_d = funEval2(T_d, delta_t, period, Q_func, RHO, dest);
   end
   if abs(F_c - F_d) < 1e-2
       search_range(cycle) = 0;
   end
end
T_ans_add = (T_serach_max+T_serach_min)/2
res = funEval2(T_ans_add, delta_t, period, Q_func, RHO, dest);
figure
plot(0:delta_t:period, P0(1:end-1))
hold on
plot([0, period], [dest, dest])
xlabel('时间(ms)')
ylabel('压力(MPa)')
title('高压油管内压力随时间变化')
plot([second, second], [0, dest])
toc
