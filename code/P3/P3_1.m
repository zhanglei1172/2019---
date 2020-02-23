%%
clc, clearvars, close all
% global data_out
data_shape = xlsread('../../data/附件1-凸轮边缘曲线.xlsx');
data_out = xlsread('../../data/附件2-针阀运动曲线.xlsx');
data_out(end,:) = [];
data = xlsread('../../data/附件3-弹性模量与压力.xlsx');
data_shape(end+1, :) = [3.14*2 data_shape(1, 2)]; % 补上2pi 便于计算
%%
omiga = 0.05545; % rad/ms
time_delay = 60.5562; % 喷嘴2的相对喷嘴1延后开启工作时间[0, 100]
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
Q_function = @(delta_P, RHO, S)C*S*sqrt(2*delta_P/RHO);
%%
delta_t = 0.1;  % 0.1ms
period_out = 100; % 100ms
% period = period_in * period_out;
period = 10e3; % 1s
% 初始化
[P0, rho, beng_P, beng_rho, sai_h]  = deal(zeros(2+period/delta_t, 1));
P0(1) = 100; % 100MPa
rho(1) = 0.850;
beng_P(1) = 0.5;
beng_rho(1) = data(data(:, 1) == beng_P(1), 3); % 高压油泵处初始密度
sai_h(1) = 0;  % 初始塞子偏移

tic
ind = uint32(1);
t = 0; % 初始时间
theta = 3.14;  % 初始角度
delta_theta = delta_t * omiga; % 每次增加的角度
% isUp = true; % 油泵是否处于上升过程
while t < period
%     gaoYaPipe_P = P0(end);
    %%%%%%%%%% 对油泵 %%%%%%%%%%
    gaoYaPipe_P = P0(ind);
    theta = theta + delta_theta;
    if theta >= 6.28
        theta = theta - 6.28;
    end
    tmp = find(data_shape(:, 1)>theta, 1);
    h = data_shape(tmp-1, 2) + 100*(theta-data_shape(tmp-1,1))*(data_shape(tmp,2)-data_shape(tmp-1,2));

    sai_h(ind+1) = h - base_h;
    delta_h = sai_h(ind+1) - sai_h(ind);
%     if abs(sai_h(ind+1) - beng_up_h) < 1e-6
%         
%     end
    if beng_P(ind) > P0(ind)
        Qin = Q_func(beng_P(ind) - gaoYaPipe_P, beng_rho(ind));
        

    else
        Qin = 0;
        if beng_P(ind) == 0 && delta_h > 0
            beng_P(ind) = beng_P(1);
            beng_rho(ind) = beng_rho(1);
        end
    end
    
    if beng_P(ind) > 0
        beng_rho(ind+1) = beng_rho(ind)*((beng_h-sai_h(ind))*beng_s - Qin*delta_t) / (beng_h-sai_h(ind+1))/beng_s; % rho = m / v
        tmp = find(beng_P(ind)<data(:, 1), 1);
        y = data(tmp-1,2) + 2*(beng_P(ind)-data(tmp-1,1))*(data(tmp,2)-data(tmp-1,2));
        beng_P(ind+1) = beng_P(ind) + (beng_rho(ind+1)-beng_rho(ind))*y/beng_rho(ind);
    end
    


    
    %%%%%%%%%% 对高压管 %%%%%%%%%%
    % 两个喷油嘴
    Qout1 = calOutQ1(mod(t, period_out), P0(ind)-0.1, rho(ind), Q_function, data_out);
    Qout2 = calOutQ1(mod(t+time_delta, period_out), P0(ind)-0.1, rho(ind), Q_function, data_out);
    Qout = Qout1 + Qout2;
%     v = delta_t * (Qin - Qout); % 增加的部分
    rho(ind+1) = rho(ind) + (Qin*beng_rho(ind)-Qout*rho(ind))*delta_t / gaoYaPipe_v;
%     P0(ind+1) = gaoYaPipe_P+ (rho(ind+1)-rho(ind))*polyval(poly, gaoYaPipe_P)/rho(ind);
    tmp = find(gaoYaPipe_P<data(:, 1), 1);
    y = data(tmp-1,2) + 2*(gaoYaPipe_P-data(tmp-1,1))*(data(tmp,2)-data(tmp-1,2));
    P0(ind+1) = gaoYaPipe_P + (rho(ind+1)-rho(ind))*y/rho(ind);

    t = t+delta_t;
    ind = ind+1;
end
figure
plot(0:delta_t:period, P0(1:end-1))
toc
hold on
plot([0, period], [100, 100])