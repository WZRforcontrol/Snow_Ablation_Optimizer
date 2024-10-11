%  MSAO: 多策略增强雪消融优化器
%
%  Developed in MATLAB R2022a
%
%  Author : Z.R.Wang
%
%  e-Mail: Email: wangzhanran@stumail.ysu.edu.cn
%
%  reference:
%  MSAO: A multi-strategy boosted snow ablation optimizer for global optimization and real-world engineering applications
%_______________________________________________________________________________________________
% 您可以简单地在一个单独的文件中定义您的成本函数，并将其句柄加载到fobj
% 您需要的初始参数如下:
%__________________________________________
% fobj = @YourCostFunction 单目标
% dim 变量个数
% Max_iteration 最大迭代次数
% SearchAgents_no 搜索代理的数目
% lb=[lb1,lb2,...,lbn] 下界
% ub=[ub1,ub2,...,ubn] 上界
% 如果所有变量的下界都相等，你可以将lb和ub定义为两个数字
%
% To run MSAO:
%______________________________________________________________________________________________

function [Best_pos,Best_score,Convergence_curve]=MSAO(N,Max_iter,lb,ub,dim,fobj)
%% hyperparameters
F0 = 0.5;%初始突变尺度因子
CR = 0.8;%交叉概率
%% 约束处理
if(max(size(ub)) == 1)
    ub = ub.*ones(1,dim);
    lb = lb.*ones(1,dim);
end
lb_matrix = repmat(lb, N, 1);
ub_matrix = repmat(ub, N, 1);

%% 良好的点集初始化策略
X=GPS_init(N,dim,ub,lb);

%% 计算第一个种群的适应度并找到最好的一个
Objective_values = zeros(1,N);
for i=1:size(X,1)
    Objective_values(i)=fobj(X(i,:));
end
[Best_score,min_ind] = min(Objective_values);
Best_pos = X(min_ind,:);

%% 构建精英池
[Obj_val_sort,idx1]=sort(Objective_values);
second_best=X(idx1(2),:);
third_best=X(idx1(3),:);
N1=floor(N*0.5);
sum1=sum(Obj_val_sort(1:N1));
half_best_mean=sum1/N1;
Elite_pool=[];
Elite_pool(1,:)=Best_pos;
Elite_pool(2,:)=second_best;
Elite_pool(3,:)=third_best;
Elite_pool(4,:)=half_best_mean;

%% other loop need parameters
Convergence_curve=[];
Convergence_curve(1) = Best_score;
index = 1:N;
Na=floor(N/2);
Nb=ceil(N/2);

%% Main loop
l=2;
while l<=Max_iter
    %% 雪融参数
    RB=randn(N,dim);%布朗随机数向量
    T=exp(-l/Max_iter);
    k=1;
    DDF=0.35*(1+(5/7)*(exp(l/Max_iter)-1)^k/(exp(1)-1)^k);
    M=DDF*T;

    %% 计算整个种群的质心位置
    X_centroid = mean(X, 1);

    %% 随机选取个体构建pop1和pop2
    index1=randperm(N,Na);
    index2=setdiff(index,index1);
    for i=1:Na
        r1=rand;
        k1=randperm(4,1);
        pa_i_old = X(index1(i),:);
        for j=1:size(X,2)
            X(index1(i),j)= Elite_pool(k1,j)+RB(index1(i),j)*(r1*(Best_pos(j)-X(index1(i),j))+(1-r1)*(X_centroid(j)-X(index1(i),j)));
        end
        % 贪心选择策略
        pa_i_val_temp = fobj(X(index1(i),:));
        if pa_i_val_temp > Objective_values(index1(i))
            X(index1(i),:) = pa_i_old;
        end
    end

    if Na<N
        Na=Na+1;
        Nb=Nb-1;
    end

    if Nb>=1
        for i=1:Nb
            r2=2*rand-1;
            pb_i_old = X(index2(i),:);
            for j=1:size(X,2)
                X(index2(i),j)= M*Best_pos(j)+RB(index2(i),j)*(r2*(Best_pos(j)-X(index2(i),j))+(1-r2)*(X_centroid(j)-X(index2(i),j)));
            end
            % 贪心选择策略
            pb_i_val_temp = fobj(X(index2(i),:));
            if pb_i_val_temp > Objective_values(index2(i))
                X(index2(i),:) = pb_i_old;
            end
        end
    end

    %% 差分进化策略
    irand = randi([1, N]);
    %% 变异操作
    F = F0*2*exp(l/(1-Max_iter));%变异尺度因子
    Mut_ref_ind = randperm(N, 3);
    V_new = X(Mut_ref_ind(1),:) + F*(X(Mut_ref_ind(2),:)-X(Mut_ref_ind(3),:));
    %% 交叉操作
    U_new = X(irand,:);
    jrand = randi([1, dim]);
    for j = 1:dim 
        if rand <= CR || j == jrand
            U_new(j) = V_new(j);
        end
    end
    %% 选择操作
    U_new_obj = fobj(U_new);
    if U_new_obj <= Objective_values(irand)
        X(irand,:) = U_new;
    end
    
    %% 检查解决方案是否超出了搜索范围，并将它们带回来
    X = max(min(X, ub_matrix), lb_matrix);

    %% 计算目标值 
    for i=1:N
        Objective_values(1,i)=fobj(X(i,:));
    end
    %% 如果有更好的解决方案，则更新目标
    [min_value,min_ind] = min(Objective_values);
    if min_value < Best_score
        Best_pos = X(min_ind,:);
        Best_score = min_value;
    else
        %% 动态种群策略
        Na = 0;
        Nb = N;
    end

    %% 动态透镜反对学习策略
    k = 1e4*(1-(l/Max_iter)^2)+1;%LOBL参数
    X(min_ind,:) = (lb + ub)./ 2 + (lb + ub) ./ (2*k) - X(min_ind,:) / k;

    %% 更新精英池
    [~,idx1]=sort(Objective_values);
    second_best=X(idx1(2),:);
    third_best=X(idx1(3),:);
    sum1=sum(Obj_val_sort(1:N1));
    half_best_mean=sum1/N1;
    Elite_pool(1,:)=Best_pos;
    Elite_pool(2,:)=second_best;
    Elite_pool(3,:)=third_best;
    Elite_pool(4,:)=half_best_mean;
    
    Convergence_curve(l)=Best_score;
    l=l+1;
end
