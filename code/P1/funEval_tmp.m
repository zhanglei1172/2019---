function f = funEval_tmp(work_time, delta_t, period, Q_func, RHO, dest)
global P0 rho gaoYaPipe_v poly weight data
% work_time = 0.291; % 单项阀每次开启时间
% delta_t = 0.01;  % 0.1ms
wait_time = 10; % 10ms
period_in = work_time + wait_time;
period_out = 100; % 100ms
% period = period_in * period_out;
% period = 1e3; % 1s
% 初始化
P_A = 160; % 160 MPa
% [P0, rho]  = deal(zeros(2+period/delta_t, 1));
P0(1) = 100; % 100MPa
rho(1) = 0.850;



ind = uint32(1);
t = 0;
while t < period
    gaoYaPipe_P = P0(ind);
    if mod(t, period_in) < work_time
        if P_A - gaoYaPipe_P < 0
            Qin = 0;
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
% f = weight(end-2e3/delta_t:end-1) * abs(dest - P0(end-2e3/delta_t:end-1));
f = mean(abs(P0(end-10000:end)-dest));
% f = -sum(log(100./P0(end-5000:end-1)));
end