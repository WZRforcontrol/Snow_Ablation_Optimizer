%%  Snow Ablation optimizer algorithm based on dynamic double population mechanism  (DP-SAO)
%
%  Developed in MATLAB R2023a
%
%  Author : Z.R.Wang
%
%  e-Mail: wangzhanran@stumail.ysu.edu.cn

%  Main paper:
%
%_______________________________________________________________________________________________
% You can simply define your cost function in a seperate file and load its handle to fobj
% The initial parameters that you need are:
%__________________________________________
% fobj = @YourCostFunction
% dim = number of your variables 决策变量个数
% Max_iteration = maximum number of iterations 最大迭代数
% SearchAgents_no = number of search agents 搜索代理的数目
% lb=[lb1,lb2,...,lbn] where lbn is the lower bound of variable n 上界
% ub=[ub1,ub2,...,ubn] where ubn is the upper bound of variable n 下界
% If all the variables have equal lower bound you can just
% define lb and ub as two single numbers

% To run DP-SAO: [Best_pos,Best_score,Convergence_curve]=DP_SAO(SearchAgents_no,Max_iteration,lb,ub,dim,fobj)
%______________________________________________________________________________________________


function [Best_pos,Best_score,Convergence_curve]=DP_SAO(N,Max_iter,lb,ub,dim,fobj)
%% 上下界避免错误
if(max(size(ub)) == 1)
    ub = ub.*ones(1,dim);
    lb = lb.*ones(1,dim);
end

%% 初始化一组随机的解 维度:𝑁×𝐷𝑖𝑚
% Initialize the set of random solutions
X=initialization_SAO(N,dim,ub,lb);
%% 相关参数运算
%% 种群所有个体适应度并计算最优个体
Objective_values = arrayfun(@(i) fobj(X(i,:)), 1:size(X,1));
[~,idx1]=sort(Objective_values);
Best_pos = X(idx1(1), :);% 最优解
Best_score = Objective_values(idx1(1));% 最优适应度
%% 精英种群𝐸𝑙𝑖𝑡𝑒(𝑡)∈[𝐺(𝑡), 𝑍𝑠𝑒𝑐𝑜𝑛𝑑 (𝑡), 𝑍𝑡ℎ𝑖𝑟𝑑 (𝑡), 𝑍𝑐(𝑡)],
Elite_pool=[];
second_best=X(idx1(2),:);
third_best=X(idx1(3),:);
N1 = floor(N * 0.5);
sum1 = sum(X(idx1(1:N1), :), 1);
half_best_mean=sum1/N1;
Elite_pool(1,:)=Best_pos;
Elite_pool(2,:)=second_best;
Elite_pool(3,:)=third_best;
Elite_pool(4,:)=half_best_mean;

%% 迭代曲线上的最优解序列
Convergence_curve=[];
Convergence_curve(1) = Best_score;

%% population
index=1:N;% 种群个体标号
%两个亚种群:𝑃𝑎和𝑃𝑏 𝑃𝑎负责开发,𝑃𝑏负责探索
Na=N/2;
Nb=N/2;

%% Main loop
% start from the second iteration since the first iteration was dedicated to calculating the fitness
% 从第二次迭代开始，因为第一次迭代专门用于计算适应度
l=2;
while l<=Max_iter
    %% 雪融参数
    RB=randn(N,dim); %Brownian random number vector 𝐵𝑀
    T=exp(-l/Max_iter);
    k=1;
    DDF=0.35*(1+(5/7)*(exp(l/Max_iter)-1)^k/(exp(1)-1)^k);
    M=DDF*T;

    %% 计算整个种群的质心位置
    % Calculate the centroid position of the entire population
    X_centroid = mean(X, 1);

    %% Select individuals randomly to construct pop1 and pop2
    % 随机选取个体构建 种群a & 种群b
    index1=randperm(N,Nb);
    index2=setdiff(index,index1);
    %% 𝑃𝑏 探索活动
    for i=1:Nb
        r1=rand;
        k1=randperm(4,1);
        for j=1:size(X,2) % in j-th dimension
            X(index1(i),j)= Elite_pool(k1,j)+RB(index1(i),j)*(r1*(Best_pos(j)-X(index1(i),j))+(1-r1)*(X_centroid(j)-X(index1(i),j)));
        end
    end

    %% 𝑃𝑎 开发活动
    for i=1:Na
        r2=2*rand-1;
        for j=1:size(X,2) % in j-th dimension
            X(index2(i),j)= M*Best_pos(j)+RB(index2(i),j)*(r2*(Best_pos(j)-X(index2(i),j))+(1-r2)*(X_centroid(j)-X(index2(i),j)));
        end
    end

    if Nb>=1
        Na=Na+rand;
        Na=round(Na);
        Nb=N-Na;
    end

    %% 检查解决方案是否超出了搜索范围，并将它们带回来
    % Check if solutions go outside the search spaceand bring them back
    for i=1:size(X,1)
        for j=1:dim
            if X(i,j)>ub(j)
                X(i,j)=ub(j);
            end
            if X(i,j)<lb(j)
                X(i,j)=lb(j);
            end
        end

        % Calculate the objective values
        Objective_values(1,i)=fobj(X(i,:));
        % Update the destination if there is a better solution
        if Objective_values(1,i)<Best_score
            Best_pos=X(i,:);
            Best_score=Objective_values(1,i);
        end
    end

    %% 种群所有个体适应度并计算最优个体
    Objective_values = arrayfun(@(i) fobj(X(i,:)), 1:size(X,1));% Calculate the objective values
    [~,idx1]=sort(Objective_values);
    Best_pos_temp = X(idx1(1), :);% 最优解
    Best_score_temp = Objective_values(idx1(1));% 最优适应度
    if Best_score_temp <= Best_score
        Best_pos = Best_pos_temp;
        Best_score = Best_score_temp;
    else
        Na=0;
        Nb=40;
    end
    %% 更新精英种群
    % Update the elite pool
    Elite_pool=[];
    second_best=X(idx1(2),:);
    third_best=X(idx1(3),:);
    N1 = floor(N * 0.5);
    sum1 = sum(X(idx1(1:N1), :), 1);
    half_best_mean=sum1/N1;
    Elite_pool(1,:)=Best_pos;
    Elite_pool(2,:)=second_best;
    Elite_pool(3,:)=third_best;
    Elite_pool(4,:)=half_best_mean;
    %% 迭代更新
    Convergence_curve(l)=Best_score;
    l=l+1;
end
