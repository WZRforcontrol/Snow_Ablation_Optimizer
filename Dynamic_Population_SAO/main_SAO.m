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
SearchAgents_no=40;
% 最大迭代次数 Maximum number of iterations   
Max_iteration=100;
% 维度(决策变量个数) Dimension,[2,10,30,50,100]
dim=10;
% 下界 lower boound
lb=-100;
% upper bound
ub=100;
% 目标函数 objective function
%
% 2020-F5，2019-F8,3 Name of the test function that can be from F1 to F10
figure
for pp = 1:9
    %% 选择目标函数
    disp(['function',num2str(pp)]);
    Function_name=pp;
    fobj = @(x) cec19_func(x',Function_name);
    %% 多次优化，消除概率数据干扰
    Max_test=10;
    for i=1:Max_test
        %% SAO优化
        [Best_pos_SAO(i,:),Best_score_SAO(i),SAO_curve_SAO(i,:)]=SAO(SearchAgents_no,Max_iteration,lb,ub,dim,fobj);
        %% DP-SAO优化
        [Best_pos_DPSAO(i,:),Best_score_DPSAO(i),SAO_curve_DPSAO(i,:)]=DP_SAO(SearchAgents_no,Max_iteration,lb,ub,dim,fobj);
    end
    %% 绘图
    subplot(3,3,pp)
    semilogy(mean(SAO_curve_SAO), 'color', 'k', 'linewidth', 2.0)
    hold on
    semilogy(mean(SAO_curve_DPSAO), 'color', 'r', 'linewidth', 2.0)
    legend('SAO', 'DP-SAO')
    f = ['Convergence curve of C2019_',num2str(pp)];
    title(f)
    xlabel('Iteration')
    ylabel('Fitness')
    axis tight
    grid off
    box on
    hold off
end
delete(MyPar)
toc
