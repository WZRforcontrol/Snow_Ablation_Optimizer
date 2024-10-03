%% Snow Ablation optimizer algorithm based on dynamic double population mechanism 
%% introduce
% Author: Z.R.Wang
% Email: wangzhanran@stumail.ysu.edu.cn
% August,5,2024 in YSU
% 动态双种群机制的雪融优化算法
% 在原算法的基础上对探索种群与开发种群数量的变化加以改进
clc;clear;close all;
rng(100)
tic
MyPar = parpool;
%% 模型超参数
% 搜索代理的数目 Number of search agents
SearchAgents_no=30;
% 最大迭代次数 Maximum number of iterations 
Max_iteration=100000;  
% 维度(决策变量个数) Dimension,[2,10,30,50,100]
dim=10; 
% 下界 lower boound
lb=-100;
% upper bound
ub=100;
% 目标函数 objective function
% 
% 2020-F5，2019-F8,3 Name of the test function that can be from F1 to F10
Function_name=8; 
fobj = @(x) cec19_func(x',Function_name);


disp('Snow ablation optimizer');
 %% SAO优化
[Best_pos_SAO,Best_score_SAO,SAO_curve_SAO]=SAO(SearchAgents_no,Max_iteration,lb,ub,dim,fobj);
%% DP-SAO优化
disp('Snow Ablation optimizer algorithm based on dynamic double population mechanism');
[Best_pos_DPSAO,Best_score_DPSAO,SAO_curve_DPSAO]=DP_SAO(SearchAgents_no,Max_iteration,lb,ub,dim,fobj);


%% 绘图
figure
dt_o = Max_iteration/10;
semilogy(SAO_curve_SAO, 'color', 'k', 'linewidth', 2.0 , 'Marker', 'o', 'MarkerIndices', 1:dt_o:length(SAO_curve_SAO))
hold on
semilogy(SAO_curve_DPSAO, 'color', 'r', 'linewidth', 2.0 , 'Marker', 'o', 'MarkerIndices', 1:dt_o:length(SAO_curve_DPSAO))
legend('SAO', 'DP-SAO')
title('Convergence curve of C2019_{8}')
xlabel('Iteration')
ylabel('Fitness')
axis tight
grid off
box on
hold off

delete(MyPar)
toc
