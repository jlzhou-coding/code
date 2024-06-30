function [OffDec,OffVel] = CompetitiveOperator2(Problem,Loser,Winner,y)
% The competitive swarm optimizer

%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Jinlong Zhou

    %% Parameter setting
    LoserDec  = Loser.decs;
    WinnerDec = Winner.decs;
    [N,D]     = size(LoserDec);
	LoserVel  = Loser.adds(zeros(N,D));


    %% Competitive swarm optimizer
    r1     = repmat(rand(N,1),1,D)/2;
    r2     = repmat(rand(N,1),1,D)/2;
    OffVel = r1.*LoserVel + r2.*(WinnerDec-LoserDec);
    OffDec = LoserDec + OffVel + r1.*(OffVel-LoserVel);
    
    %% Polynomial mutation
    Lower  = repmat(Problem.lower,N,1);
    Upper  = repmat(Problem.upper,N,1);
    disM   = 20;
    Site   = rand(N,D) < 1/D;
    mu     = rand(N,D);
    temp   = Site & mu<=0.5;
    OffDec       = max(min(OffDec,Upper),Lower);
    OffDec(temp) = OffDec(temp)+(Upper(temp)-Lower(temp)).*((2.*mu(temp)+(1-2.*mu(temp)).*...
                   (1-(OffDec(temp)-Lower(temp))./(Upper(temp)-Lower(temp))).^(disM+1)).^(1/(disM+1))-1);
    temp  = Site & mu>0.5; 
    OffDec(temp) = OffDec(temp)+(Upper(temp)-Lower(temp)).*(1-(2.*(1-mu(temp))+2.*(mu(temp)-0.5).*...
                   (1-(Upper(temp)-OffDec(temp))./(Upper(temp)-Lower(temp))).^(disM+1)).^(1/(disM+1)));
end