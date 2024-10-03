clc;clear;close all;
load('grouped_data.mat');
load('wwwwwwwwwwqx.mat')
n = height(grouped_data);
%% 
x = optimvar("x",1,n,"LowerBound",0);
x_is = optimvar("x_is",1,n,"Type","integer","LowerBound",0,"UpperBound",1);
k_u_d = optimvar("k_u_d",2,n,"LowerBound",0);

initialPoint.x = zeros(size(x));
initialPoint.x_is = zeros(size(x_is));
initialPoint.k_u_d = zeros(size(k_u_d));

prob1 = optimproblem;

prob1.Objective = sum(k_u_d,'all');

prob1.Constraints.constraint1 = sum(x_is) <= 33;
prob1.Constraints.constraint2 = sum(x_is) >= 27;
prob1.Constraints.constraint3 = x.*x_is <= x;
prob1.Constraints.constraint4 = x >= x_is*2.5;
prob1.Constraints.constraint5 = x -k_u_d(1,:)+ k_u_d(2,:) == wqxxxxxxxxxxxxxxxx;

options = optimoptions("ga","Display","iter","PlotFcn",@WZR_gaplotbestf);

show(prob1);

[sol1,fval1,~] = solve(prob1,initialPoint,"Solver","ga","Options",options);
fig = gcf;
saveas(fig, 'ga_iteration_plot1.png');

prob2 = optimproblem("ObjectiveSense","Maximize");
prob2.Constraints.constraint1 = sum(x_is) <= 33;
prob2.Constraints.constraint2 = sum(x_is) >= 27;
prob2.Constraints.constraint3 = x.*x_is <= x;
prob2.Constraints.constraint4 = x >= x_is*2.5;
prob2.Constraints.constraint5 = x -k_u_d(1,:)+ k_u_d(2,:) == wqxxxxxxxxxxxxxxxx;
prob2.Constraints.constraint6 = sum(k_u_d,'all') <= fval1;

prob2.Objective = x*(grouped_data.danjia-grouped_data.pifajia);

options = optimoptions("ga","Display","iter","PlotFcn",@WZR_gaplotbestf);

show(prob2);

[sol2,fval2,~] = solve(prob2,sol1,"Solver","ga","Options",options);
fig = gcf;
saveas(fig, 'ga_iteration_plot2.png');

for i = 1:n
    if sol2.x_is(i) == 1
        aaqaaaaa = cellstr(grouped_data.danpinname(i));
        aaqaaaaa = [aaqaaaaa  ' '];
        disp([aaqaaaaa{1},num2str(sol2.x(i)),'kg'])
    end
end
disp(['家人们，咱们可以搞',num2str(fval2),'元'])


