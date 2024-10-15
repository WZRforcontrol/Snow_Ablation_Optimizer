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
% To run MSAO: [Best_pos, Best_score, Convergence_curve] = MSAO(N, Max_iter, lb, ub, dim, fobj)
%______________________________________________________________________________________________

function [Best_pos, Best_score, Convergence_curve] = MSAO(N, Max_iter, lb, ub, dim, fobj)
    %% 超参数设置
    F0 = 0.5; % 初始变异因子
    CR = 0.8; % 交叉概率

    %% 约束处理
    if numel(ub) == 1
        ub = ub * ones(1, dim);
        lb = lb * ones(1, dim);
    end

    %% 良好的点集初始化策略
    X = GPS_init(N, dim, ub, lb);

    %% 计算初始种群的适应度并找到最佳个体
    Objective_values = zeros(1, N);
    for i = 1:N
        Objective_values(i) = fobj(X(i, :));
    end
    [Best_score, min_ind] = min(Objective_values);
    Best_pos = X(min_ind, :);

    %% 构建精英池
    [~, idx1] = sort(Objective_values);
    second_best = X(idx1(2), :);
    third_best = X(idx1(3), :);
    N1 = floor(N * 0.5);
    half_best_mean = mean(X(idx1(1:N1), :), 1);
    Elite_pool = [];
    Elite_pool(1, :) = Best_pos;
    Elite_pool(2, :) = second_best;
    Elite_pool(3, :) = third_best;
    Elite_pool(4, :) = half_best_mean;

    %% 初始化其他参数
    Convergence_curve = zeros(1, Max_iter);
    Convergence_curve(1) = Best_score;
    index = 1:N;
    Na = floor(N / 2);
    Nb = ceil(N / 2);

    %% 主循环
    l = 2;
    while l <= Max_iter
        %% 雪融参数计算
        RB = randn(N, dim); % 布朗运动随机数向量
        T = exp(-l / Max_iter);
        DDF = 0.35 + 0.25 * (exp(l / Max_iter) - 1) / (exp(1) - 1);
        M = DDF * T;

        %% 计算种群质心
        X_centroid = mean(X, 1);

        %% 构建子种群Pa和Pb
        index1 = randperm(N, Na);
        index2 = setdiff(index, index1);

        %% 探索阶段（子种群Pa）
        for i = 1:Na
            r1 = rand;
            k1 = randi([1, 4], 1);
            pa_i_old = X(index1(i), :);
            X_new = Elite_pool(k1, :) + RB(index1(i), :) .* ...
                    (r1 * (Best_pos - pa_i_old) + (1 - r1) * (X_centroid - pa_i_old));
            % 边界检查
            X_new = max(min(X_new, ub), lb);
            % 计算新适应度值
            pa_i_val_temp = fobj(X_new);
            % 贪心选择策略
            if pa_i_val_temp < Objective_values(index1(i))
                X(index1(i), :) = X_new;
                Objective_values(index1(i)) = pa_i_val_temp;
            end
            % 否则保持原来的值
        end

        %% 动态调整子种群大小
        if Na < N
            Na = Na + 1;
            Nb = Nb - 1;
            if Nb < 0
                Nb = 0;
            end
        end

        %% 开发阶段（子种群Pb）
        if Nb >= 1
            for i = 1:Nb
                r2 = 2 * rand - 1;
                pb_i_old = X(index2(i), :);
                X_new = M * Best_pos + RB(index2(i), :) .* ...
                        (r2 * (Best_pos - pb_i_old) + (1 - r2) * (X_centroid - pb_i_old));
                % 边界检查
                X_new = max(min(X_new, ub), lb);
                % 计算新适应度值
                pb_i_val_temp = fobj(X_new);
                % 贪心选择策略
                if pb_i_val_temp < Objective_values(index2(i))
                    X(index2(i), :) = X_new;
                    Objective_values(index2(i)) = pb_i_val_temp;
                end
                % 否则保持原来的值
            end
        end

        %% 差分进化策略
        for irand = 1:N
            %% 变异操作
            F = F0 * 2 * exp(-l / Max_iter); % 变异尺度因子
            Mut_ref_ind = randperm(N, 3);
            while any(Mut_ref_ind == irand)
                Mut_ref_ind = randperm(N, 3);
            end
            V_new = X(Mut_ref_ind(1), :) + F * (X(Mut_ref_ind(2), :) - X(Mut_ref_ind(3), :));
            % 边界检查
            V_new = max(min(V_new, ub), lb);
            %% 交叉操作
            U_new = X(irand, :);
            jrand = randi([1, dim]);
            for j = 1:dim
                if rand <= CR || j == jrand
                    U_new(j) = V_new(j);
                end
            end
            % 边界检查
            U_new = max(min(U_new, ub), lb);
            %% 选择操作
            U_new_obj = fobj(U_new);
            if U_new_obj < Objective_values(irand)
                X(irand, :) = U_new;
                Objective_values(irand) = U_new_obj;
            end
            % 否则保持原来的值
        end

        %% 更新全局最优解
        [min_value, min_ind_new] = min(Objective_values);
        if min_value <= Best_score
            Best_pos = X(min_ind_new, :);
            Best_score = min_value;
            min_ind = min_ind_new;
        else
            % 找到Best_pos在X中的索引
            min_ind = find(ismember(X, Best_pos, 'rows'), 1);
            if isempty(min_ind)
                % 如果Best_pos不在X中，则跳过更新
                min_ind = 1; % 默认设置为1，防止后续报错
            end
        end

        %% 动态透镜反对学习策略
        k_LOBL = 1e4 * (1 - (l / Max_iter)^2) + 1; % LOBL参数
        x_temp = (lb + ub) / 2 + ((lb + ub) / (2 * k_LOBL) - Best_pos / k_LOBL);
        % 边界检查
        x_temp = max(min(x_temp, ub), lb);
        fit_temp = fobj(x_temp);
        if fit_temp < Best_score
            Best_pos = x_temp;
            Best_score = fit_temp;
            % 不更新X和Objective_values
        end

        %% 更新精英池
        [~, idx1] = sort(Objective_values);
        second_best = X(idx1(2), :);
        third_best = X(idx1(3), :);
        half_best_mean = mean(X(idx1(1:N1), :), 1);
        Elite_pool = [];
        Elite_pool(1, :) = Best_pos;
        Elite_pool(2, :) = second_best;
        Elite_pool(3, :) = third_best;
        Elite_pool(4, :) = half_best_mean;

        %% 记录收敛曲线
        Convergence_curve(l) = Best_score;
        l = l + 1;
    end
end