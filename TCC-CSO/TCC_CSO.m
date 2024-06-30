classdef TCC_CSO < ALGORITHM
% <multi> <real> <large/none> <constrained>
% competitive swarm algorithm based on tri-competitive criterion

%------------------------------- Reference --------------------------------
% J. Zhou, Y. Zhang, and Y. Fan, “A competitive swarm algorithm based on tri-competitive criterion for
% constrained multiobjective optimization,” in Genetic and Evolutionary Computation Conference
% (GECCO’24 Companion), July 14–18, 2024, Melbourne, VIC, Australia. ACM, New York, NY, USA, 4 pages.  
%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Jinlong Zhou

    methods
        function main(Algorithm,Problem)
            %% Generate random population
            Population = Problem.Initialization();
            CV         = sum(max(0,Population.cons),2);
            CVmax      = max(CV);
            epsilon_0  = CVmax;
            epsilon    = epsilon_0;
            Competitive_Pop = UpdateP1(Population,Problem.N,epsilon);
            [Cooperative_Pop,Cooperative_Pop_Fitness] = UpdateP2(Population,Problem.N);
            Population = ArchiveUpdate(Population,Problem.N);
            Tc    = 0.9 * ceil(Problem.maxFE/Problem.N);
            cp    = 2;
            alpha = 0.95;
            tao   = 0.05;
            y     = 10;
            G     = Problem.maxFE/Problem.N;
            
            %% Optimization
            while Algorithm.NotTerminated(Population)
                gen       = ceil(Problem.FE/Problem.N);
                CV        = sum(max(0,Competitive_Pop.cons),2);
                CV_max    = max(CV);
                CVmax     = max([CV_max,CVmax]);
                epsilon_0 = CVmax;
                rf = sum(CV <= 1e-6) / length(Competitive_Pop);
                epsilon = update_epsilon(tao,epsilon,epsilon_0,rf,alpha,gen,Tc,cp);
                
                Competitive_Pop_Fitness = CalFitness(Competitive_Pop.objs,Competitive_Pop.cons,epsilon);
                
                if length(Competitive_Pop) >= 3
                    Rank = randperm(length(Competitive_Pop),floor(length(Competitive_Pop)/3)*3);
                else
                    Rank = [1,1,1];
                end
                Loser1  = Rank(1:end/3);
                Loser2  = Rank(end/3+1:end/3*2);
                Winner = Rank(end/3*2+1:end);
                 %% group winner and loser
                Change = Competitive_Pop_Fitness(Loser2) <= Competitive_Pop_Fitness(Winner);
                Temp   = Winner(Change);
                Winner(Change) = Loser2(Change);
                Loser2(Change)  = Temp;%group winner and loser

                ChangeT = Competitive_Pop_Fitness(Loser1) <= Competitive_Pop_Fitness(Loser2);
                Temp   = Loser2(ChangeT);
                Loser2(ChangeT) = Loser1(ChangeT);
                Loser1(ChangeT)  = Temp;%group winner and loser
                
                [Loser2Dec,Loser2Vel] = CompetitiveOperator1(Problem,Competitive_Pop(Loser2),Competitive_Pop(Winner),y);
                Loser2_Pop = Problem.Evaluation(Loser2Dec,Loser2Vel);
         
                [Loser1Dec,Loser1Vel] =  CompetitiveOperator2(Problem,Competitive_Pop(Loser1),Loser2_Pop,y);%%%
                Loser1_Pop = Problem.Evaluation(Loser1Dec,Loser1Vel);
                Offspring1 = [Loser1_Pop,Loser2_Pop,Competitive_Pop(Winner)];
                
                LearningPool = TournamentSelection(2,Problem.N,Cooperative_Pop_Fitness);
                Offspring2   = CooperativeOperator(Problem,Cooperative_Pop(LearningPool));
                
                Offspring = [Offspring1,Offspring2];
                
                Population = ArchiveUpdate([Population,Offspring],Problem.N);
                Competitive_Pop = UpdateP1([Competitive_Pop,Offspring],Problem.N,epsilon);
                gen = ceil(Problem.FE/Problem.N);
                y   = ((Problem.M)^2+1)-((Problem.M)^2)./(1+exp(-10*(gen/G-0.5)));
                [Cooperative_Pop,Cooperative_Pop_Fitness] = UpdateP2([Offspring,Cooperative_Pop],Problem.N);
            end
        end
    end
end

function epsilon = update_epsilon(tao,epsilon_k,epsilon_0,rf,alpha,gen,Tc,cp)
% Update the value of epsilon

    if gen > Tc
        epsilon = 0;
    else
        if rf < alpha
            epsilon = (1 - tao) * epsilon_k;
        else
            epsilon = epsilon_0 * ((1 - (gen / Tc)) ^ cp);
        end
    end
end