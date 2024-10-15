
% 该函数通过佳点集初始化策略初始化第一个搜索代理种群
function X = GPS_init(N, dim, ub, lb)
    % GPS_init - Good Point Set initialization
    % N - 种群规模
    % dim - 变量维度
    % ub - 上界（向量）
    % lb - 下界（向量）
    
    % 找到满足条件的最小素数p
    p_candidate = 2 * dim + 3;
    while ~isprime(p_candidate)
        p_candidate = p_candidate + 1;
    end
    p = p_candidate;
    
    % 计算r_j值
    r = zeros(1, dim);
    for j = 1:dim
        r(j) = mod(2 * cos(2 * pi * j / p), 1);
    end
    
    % 构造N个点
    P = zeros(N, dim);
    for i = 1:N
        for j = 1:dim
            P(i, j) = mod(r(j) * i, 1);
        end
    end
    
    % 映射到可行域
    X = lb + P .* (ub - lb);
end



