function [ Q ] = calOutQ( time_in_period, delta_p, rho, Q_function )
%CALOUTQ 此处显示有关此函数的摘要
%   此处显示详细说明
global data_out
if time_in_period > 2.45
    Q = 0;
    return
elseif time_in_period > 2.000001
    tmp = find(time_in_period<=data_out(:, 4), 1);
    if tmp == 1
        h = 1.9971;
    else
        h = data_out(tmp-1,5) + 100*(time_in_period-data_out(tmp-1,4))*(data_out(tmp,5)-data_out(tmp-1,5));

    end
    
elseif time_in_period < 0.45
    tmp = find(time_in_period>=data_out(:, 1), 1, 'last')+1;
    if tmp > 45
        h = 1.9975;
    else
        h = data_out(tmp-1,2) + 100*(time_in_period-data_out(tmp-1,1))*(data_out(tmp,2)-data_out(tmp-1,2));
    end
    
else
    h = 2;
    
end
if h == 0
    Q = 0;
    return
end
s = min([pi*( (tand(9)*(1.25/tand(9) + h))^2 - 1.25^2), 1.5394]); % tan9 = 0.1584
Q = Q_function(delta_p, rho, s);
end

