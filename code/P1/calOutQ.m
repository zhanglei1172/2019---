function [ Q ] = calOutQ( time_in_period )
%CALOUTQ 此处显示有关此函数的摘要
%   此处显示详细说明
if time_in_period > 2.4
    Q = 0;
elseif time_in_period > 2.2
    Q = -100* time_in_period + 240;
elseif time_in_period < 0.2
    Q = 100 * time_in_period;
else
    Q = 20;
end
end

