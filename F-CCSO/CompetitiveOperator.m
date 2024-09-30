function [Off_P,Off_V] = CompetitiveOperator(Problem,Population,fit)
% The competitive swarm optimizer
P_Dec = Population.decs;
[N,D] = size(P_Dec);
V     = Population.adds(zeros(N,D));
Off_P = zeros(N,D);
Off_V = zeros(N,D);

%% Learning
[~,rank]  =sort(fit);
r=Problem.FE/Problem.maxFE;
X=r;%line
Winner_N=round(X*20);

LeaderSet = rank(1:max(Winner_N,2));
for i = 1 : N% update each particle
    winner = LeaderSet(randperm(length(LeaderSet),2));%
    mask   = (fit(winner(1)) > fit(winner(2)));
    winner = ~mask.*winner(1) + mask.*winner(2);% the smaller fitness is the winner
    % Learning:  i learn from the winner
    r1 = rand(1,D);
    r2 = rand(1,D);
    Off_V(i,:) = r1.*V(i,:) + r2.*(P_Dec(winner,:)-P_Dec(i,:));%*y;%update V
    Off_P(i,:) = P_Dec(i,:) + Off_V(i,:) + r1.*(Off_V(i,:)-V(i,:));%*(-1)^n;
end

%% Polynomial mutation
proM=1;
disM=20;
Lower = repmat(Problem.lower,N,1);
Upper = repmat(Problem.upper,N,1);
Site  = rand(N,D) < proM/D;
mu    = rand(N,D);
temp  = Site & mu<=0.5;
Off_P = min(max(Off_P,Lower),Upper);
Off_P(temp) = Off_P(temp)+(Upper(temp)-Lower(temp)).*((2.*mu(temp)+(1-2.*mu(temp)).*...
    (1-(Off_P(temp)-Lower(temp))./(Upper(temp)-Lower(temp))).^(disM+1)).^(1/(disM+1))-1);
temp = Site & mu>0.5;
Off_P(temp) = Off_P(temp)+(Upper(temp)-Lower(temp)).*(1-(2.*(1-mu(temp))+2.*(mu(temp)-0.5).*...
    (1-(Upper(temp)-Off_P(temp))./(Upper(temp)-Lower(temp))).^(disM+1)).^(1/(disM+1)));
end