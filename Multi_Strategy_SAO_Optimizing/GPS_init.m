
% �ú���ͨ���ѵ㼯��ʼ�����Գ�ʼ����һ������������Ⱥ
function X = GPS_init(N, dim, ub, lb)
    % GPS_init - Good Point Set initialization
    % N - ��Ⱥ��ģ
    % dim - ����ά��
    % ub - �Ͻ磨������
    % lb - �½磨������
    
    % �ҵ�������������С����p
    p_candidate = 2 * dim + 3;
    while ~isprime(p_candidate)
        p_candidate = p_candidate + 1;
    end
    p = p_candidate;
    
    % ����r_jֵ
    r = zeros(1, dim);
    for j = 1:dim
        r(j) = mod(2 * cos(2 * pi * j / p), 1);
    end
    
    % ����N����
    P = zeros(N, dim);
    for i = 1:N
        for j = 1:dim
            P(i, j) = mod(r(j) * i, 1);
        end
    end
    
    % ӳ�䵽������
    X = lb + P .* (ub - lb);
end



