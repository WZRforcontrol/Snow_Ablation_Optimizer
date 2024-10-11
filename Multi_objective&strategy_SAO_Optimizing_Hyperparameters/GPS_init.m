
% �ú���ͨ���ѵ㼯��ʼ�����Գ�ʼ����һ������������Ⱥ
function Positions=GPS_init(SearchAgents_no,dim,ub,lb)
% SearchAgents_no����Ⱥ����
% dim��ά��
% ub lb��ȡֵ��Χ
Boundary_no= size(ub,2); % numnber of boundaries

% ������б����ı߽綼��ȣ������û�Ϊub��l������һ������
if Boundary_no==1 && dim ~= 1
    lb = lb * ones(1,dim);
    ub = ub * ones(1,dim);
%     Positions=rand(SearchAgents_no,dim).*(ub-lb)+lb;
end
p = zeros(SearchAgents_no,dim);
prime_number_min = dim*2 +3;
% �ҵ�(prime_number_min-3)/2>=dim����С����prime_number_min
while 1
    if isprime(prime_number_min)==1
        break;
    else
       prime_number_min = prime_number_min + 1;
    end
end

for i = 1:SearchAgents_no
    for j = 1:dim
        r = mod(2*cos(2*pi*j/prime_number_min)*i,1);% ��Ӧά�ȵ�r
%         r = mod(exp(j)*i,1);
        p(i,j) = lb(j)+r*(ub(j)-lb(j));
    end
end
Positions = p;
end


