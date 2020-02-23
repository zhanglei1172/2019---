%% 灵敏度分析：对最优omiga波动1%造成的压强影响
clc, clearvars, close all
global data_out data_shape data sai_h beng_rho beng_P P0 base_h beng_h beng_s period_out
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
delta_t = 0.1;  % 0.1ms
period_out = 100; % 100ms
% period = period_in * period_out;
period = 10e3; % 1s
% 初始化
global rho weight
[P0, rho, beng_P, beng_rho, sai_h]  = deal(zeros(2+period/delta_t, 1));
P0(1) = 100; % 100MPa
rho(1) = 0.850;
beng_P(1) = 0.5;
beng_rho(1) = data(data(:, 1) == beng_P(1), 3); % 高压油泵处初始密度
sai_h(1) = 0;  % 初始塞子偏移
omiga = 0.02767; % rad/ms
tic

subplot(311)
funEval(omiga*0.999, delta_t, period, Q_func);
plot(0:delta_t:period, P0(1:end-1))
hold on
plot([0, period], [100, 100])
mean(abs(P0(end-10000:end-1)-100))
toc

funEval(omiga, delta_t, period, Q_func);
subplot(312)
plot(0:delta_t:period, P0(1:end-1))
hold on
plot([0, period], [100, 100])
mean(abs(P0(end-10000:end-1)-100))
toc


funEval(omiga*1.001, delta_t, period, Q_func);
subplot(313)
plot(0:delta_t:period, P0(1:end-1))
hold on
plot([0, period], [100, 100])
mean(abs(P0(end-10000:end-1)-100))
toc
