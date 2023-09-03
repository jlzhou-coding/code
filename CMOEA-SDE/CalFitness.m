function Fitness = CalFitness(PopObj,CV,r)
% Calculate the fitness of each solution

%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------
    N = size(PopObj,1);
    %% Detect the dominance relation between each two solutions
    Dominate = false(N);
    for i = 1 : N-1
        for j = i+1 : N
            k = any(PopObj(i,:)<PopObj(j,:)) - any(PopObj(i,:)>PopObj(j,:));
            if k == 1 && 1*CV(i)<=CV(j)
                Dominate(i,j) = true;
            elseif k == -1 && CV(i)>=CV(j)*1
                Dominate(j,i) = true;
            end
        end
    end

    %% Calculate S(i)：i支配的个体数
    S = sum(Dominate,2);

    %% Calculate R(i):支配i个体集合S，所支配的个体数总和
    R = zeros(1,N);
    for i = 1 : N
        R(i) = sum(S(Dominate(:,i)));
    end

    sigmoid=1.0./(1+exp(-25*(r-0.25)));
    %% Calculate the shifted distance between each two solutions

    Distance = inf(N);
    for i = 1 : N
        SPopObj0= PopObj;
        for k=1:N
            if CV(i)>CV(k)
                SPopObj(k,:)  =SPopObj0(k,:)+(PopObj(i,:)-SPopObj0(k,:))*sigmoid;
            else
                SPopObj(k,:)=SPopObj0(k,:);
            end
        end

        for j = [1:i-1,i+1:N]
            Distance(i,j) = norm(PopObj(i,:)-SPopObj(j,:));
        end
    end

    %% Calculate D(i)
    Distance = sort(Distance,2);
    D = 1./(Distance(:,floor(sqrt(N)))+2);

    %% Calculate the fitnesses
    Fitness = R + D';
end