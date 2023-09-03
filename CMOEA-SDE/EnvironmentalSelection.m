function [Population,Fitness] = EnvironmentalSelection(Population,N,r)
% The environmental selection of SPEA2+CSDE

%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------
    %% Calculate the CV of each solution
    CV=sum(max(Population.cons,0),2);
    %% Calculate the fitness of each solution
    Fitness = CalFitness(Population.objs,CV,r);
  
    %% Environmental selection
    Next = Fitness < 1;
    if sum(Next) < N
        [~,Rank] = sort(Fitness);
        Next(Rank(1:N)) = true;
    elseif sum(Next) > N
        Del  = Truncation(Population(Next).objs,sum(Next)-N,CV(Next),r);
        Temp = find(Next);
        Next(Temp(Del)) = false;
    end
    % Population for next generation
    Population = Population(Next);
    Fitness    = Fitness(Next);
end

function Del = Truncation(PopObj,K,CV,r)
% Select part of the solutions by truncation
    N = size(PopObj,1);
    %% Calculate the shifted distance between each two solutions
    Distance = inf(N);
    sigmoid=1.0./(1+exp(-25*(r-0.25)));
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

    %% Truncation
    Del = false(1,N);
    while sum(Del) < K
        Remain   = find(~Del);
        Temp     = sort(Distance(Remain,Remain),2);
        [~,Rank] = sortrows(Temp);
        Del(Remain(Rank(1))) = true;
    end
end