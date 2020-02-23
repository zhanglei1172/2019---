function f = funEval(omiga, delta_t, period, Q_func, time_delta, data_out)
time_delta = 100 - time_delta;
global P0 rho gaoYaPipe_v data data_shape sai_h beng_rho beng_P base_h beng_h beng_s period_out Q_function
beng_P  = zeros(2+period/delta_t, 1);
beng_P(1) = 0.5;
ind = uint32(1);
t = 0; % 初始时间
theta = 3.14;  % 初始角度
delta_theta = delta_t * omiga; % 每次增加的角度
% isUp = true; % 油泵是否处于上升过程
% P0(1) = 100; % 100MPa
% rho(1) = 0.850;
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
    
    
    Qout1 = calOutQ1(mod(t, period_out), P0(ind)-0.1, rho(ind), Q_function, data_out);
    Qout2 = calOutQ1(mod(t+time_delta, period_out), P0(ind)-0.1, rho(ind), Q_function, data_out);
    Qout = Qout1 + Qout2;
%     Qout = calOutQ(mod(t, period_out), P0(ind)-0.1, rho(ind), Q_function);
%     v = delta_t * (Qin - Qout); % 增加的部分
    rho(ind+1) = rho(ind) + (Qin*beng_rho(ind)-Qout*rho(ind))*delta_t / gaoYaPipe_v;
%     P0(ind+1) = gaoYaPipe_P+ (rho(ind+1)-rho(ind))*polyval(poly, gaoYaPipe_P)/rho(ind);
    tmp = find(gaoYaPipe_P<data(:, 1), 1);
    y = data(tmp-1,2) + 2*(gaoYaPipe_P-data(tmp-1,1))*(data(tmp,2)-data(tmp-1,2));
    P0(ind+1) = gaoYaPipe_P + (rho(ind+1)-rho(ind))*y/rho(ind);

    t = t+delta_t;
    ind = ind+1;
end
% f = weight(end-2e3/delta_t:end-1) * abs(100 - P0(end-2e3/delta_t:end-1));
f = mean(abs(P0(end-10000:end-1)-100));
end