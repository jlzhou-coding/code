function [Population,FrontNo]= EnvironmentalSelection(Population,N)
%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

%% Non-dominated sorting
[FrontNo,MaxFNo] = NDSort(Population.objs,N);

%% Association operation
Next   = [find(FrontNo<MaxFNo),find(FrontNo==MaxFNo)];
Choose = Association(Population(FrontNo<MaxFNo).objs,Population(FrontNo==MaxFNo).objs,N);
Next   = Next(Choose);
% Population for next generation
Population = Population(Next);
FrontNo    = FrontNo(Next);
end

function Choose = Association(PopObj1,PopObj2,N)
% Association operation in the algorithm

[N1,~] = size(PopObj1);
[N2,M] = size(PopObj2);
PopObj = [PopObj1;PopObj2];

%% Normalization
Zmin   = min(PopObj,[],1);
Zmax   = max(PopObj,[],1);
PopObj = (PopObj-repmat(Zmin,size(PopObj,1),1))./repmat(Zmax-Zmin,size(PopObj,1),1);

%% Calculate the fitness value of each solution
fit = sqrt(sum(PopObj.^2,2));
angle = acos(1-pdist2(PopObj,PopObj,'cosine'));
%% Niching
Choose = [true(1,N1),false(1,N2)];


Del=false(1,N1+N2);
angle = acos(1-pdist2(PopObj,PopObj,'cosine'));
angle(logical(eye(length(angle)))) = inf;
sa=sort(angle,2);
r  = median(sa(:,min(M,size(sa,2))));
R  = min(angle./r,1);


DI=1-prod(R,2);

temp=pdist2(PopObj,eye(M),'cosine');
[~,extreme]= min(temp,[],1);
DI(extreme)=-1;
fit(extreme)=0;
% Delete solution one by one
while sum(Del) < N1+N2-N
    Remain=find(~Del & ~Choose);
    DI2=DI(Remain);
    [~,worst]  = max(DI2);
    individual=Remain(worst);
    Remain(worst)=[];
    
    [~,I]=find(R(individual,:)<1) ;
    I=intersect(I,Remain);
    if ~isempty(I)
        [theta,r]=min(angle(individual,I));
        if fit(I(r))>fit(individual)
            individual2=I(r);
        else
            individual2=individual;
        end
    else
        individual2=individual;
    end
    
    Del(individual2)=true;
    
    [~,I_niche]=find(R(individual2,:)<1) ;
    
    R(individual2,:) =1;
    R(:,individual2)=1;
    
    
    DI(I_niche)=1-prod(R(I_niche,:),2);
    DI(individual2)=-1;
    DI(extreme)=-1;
end
Choose=~Del;
end