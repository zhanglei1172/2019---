clc, clearvars, close all
global data_out
data_shape = xlsread('../../data/附件1-凸轮边缘曲线.xlsx');
data_out = xlsread('../../data/附件2-针阀运动曲线.xlsx');
data_out(end,:) = [];
data = xlsread('../../data/附件3-弹性模量与压力.xlsx');

%% 凸轮曲线
% reglm(data_shape(:, 2),cos(data_shape(:, 1)))
func = @(p, x)p(2)+p(1)*cos(x);
cof = lsqcurvefit(func, [1,1],data_shape(:, 1), data_shape(:, 2))

plot(data_shape(:, 1), data_shape(:, 2), 'bo')
hold on
plot(data_shape(:, 1), func(cof, data_shape(:, 1)), 'r-', 'LineWidth', 2)

xlabel('\theta')
ylabel('H')
%% 针阀升呈
% func = @(a, x)a(1)./(1+(a(1)/a(2)-1)*exp(-a(3).*x));
% a=lsqcurvefit(func,[1, 1, 1],data_out(:, 1),data_out(:, 2));
% plot(data_out(:, 1),func(a, data_out(:, 1)))
% hold on
% plot(data_out(:, 1),data_out(:, 2))
% fun = @(p, x)p(1)*exp(p(2)./(x+p(3)));
% p0 = [10, 1, 1];
% opt = statset;
% opt.Robust = 'on';
% model = NonLinearModel.fit(data_out(:, 1),data_out(:, 2),fun,p0,'Options',opt)
%% 弹性模量
figure
poly = polyfit(data(:, 1), data(:, 2), 3);
plot(data(:, 1), data(:, 2), '*')
hold on % 压力 -> 弹性模量
plot(0:0.1:200, polyval(poly, 0:0.1:200), 'LineWidth', 2)
xlabel('E')
ylabel('P')
disp([num2str(poly(1)), 'x^3','+',num2str(poly(2)), 'x^2','+', num2str(poly(3)), 'x','+', num2str(poly(4))])