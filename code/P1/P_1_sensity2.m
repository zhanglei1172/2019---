%% 灵敏度分析：在最优开启时间下对高压油泵恒压160MPa上下波动1Mpa的影响
clc, clearvars, close all;
data = xlsread('../../data/附件3-弹性模量与压力.xlsx');
poly = polyfit(data(:, 1), data(:, 2), 4);
plot(data(:, 1), data(:, 2), '*')
hold on % 压力 -> 弹性模量
plot(0:0.1:200, polyval(poly, 0:0.1:200), 'LineWidth', 2)
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
RHO = data(data(:, 1) == 160, 3); % 高压油泵处密度
%% 设置参数
gaoYaPipe_length = 500; % 500 mm长
gaoYaPipe_diameter = 10;
gaoYaPipe_v = gaoYaPipe_length * (pi*(gaoYaPipe_diameter/2)^2); % 体积
A_diameter = 1.4;
% P_A = 160; % 160 MPa
% 单向阀
wait_time = 10; % 10ms

% 喷油器
dest = 100;
% 公式参数
C = 0.85;
A = pi*(A_diameter/2)^2;
Q_func = @(delta_P)C*A*sqrt(2*delta_P/RHO);
%%
% tic
% profile on
work_times = 0.2877; % 单项阀每次开启时间
work_time = work_times;
delta_t = 0.1;  % 0.1ms
period_in = work_time + wait_time;
period_out = 100; % 100ms
% period = period_in * period_out;
period = 10e3; % 1s
% 初始化
global P0
[P0, rho]  = deal(zeros(2+period/delta_t, 1));
P0(1) = 100; % 100MPa
rho(1) = 0.850;

% P_A = 160; % 160 MPa
cnt = 0;
for P_A = 159:1:161
ind = uint32(1);
t = 0;
while t < period
%     gaoYaPipe_P = P0(end);
    gaoYaPipe_P = P0(ind);
    if mod(t, period_in) < work_time
        if P_A - gaoYaPipe_P < 0
            Qin = 0;
%             return;
        else
            Qin = Q_func(P_A - gaoYaPipe_P); % P_A > P 压力差
        end
    else
        Qin = 0;
    end
    
    Qout = calOutQ(mod(t, period_out)); 
%     v = delta_t * (Qin - Qout); % 增加的部分
    rho(ind+1) = rho(ind) + (Qin*RHO-Qout*rho(ind))*delta_t / gaoYaPipe_v;
%     P0(ind+1) = gaoYaPipe_P+ (rho(ind+1)-rho(ind))*polyval(poly, gaoYaPipe_P)/rho(ind);
    tmp = find(gaoYaPipe_P<data(:, 1), 1);
    y = data(tmp-1,2) + 2*(gaoYaPipe_P-data(tmp-1,1))*(data(tmp,2)-data(tmp-1,2));
    P0(ind+1) = gaoYaPipe_P + (rho(ind+1)-rho(ind))*y/rho(ind);

    t = t+delta_t;
    ind = ind+1;
end
subplot(311+cnt)
mean(abs(P0(end-5000:end)-dest))
plot(0:delta_t:period, P0(1:end-1))
title(num2str(P_A), 'FontSize',10)
hold on
plot([0, period], [100, 100])
xlabel('时间(ms)')
ylabel('压力(MPa)')
cnt = cnt+1;
end


% toc
% 
% funEval(work_time*0.99, delta_t, period, Q_func, RHO, dest);
% subplot(312)
% plot(0:delta_t:period, P0(1:end-1))
% hold on
% plot([0, period], [100, 100])
% toc
% 
% 
% funEval(work_time*1.01, delta_t, period, Q_func, RHO, dest);
% subplot(313)
% plot(0:delta_t:period, P0(1:end-1))
% hold on
% plot([0, period], [100, 100])
% toc
