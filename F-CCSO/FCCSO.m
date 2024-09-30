classdef FCCSO < ALGORITHM
% <multi> <real> <large/none> <constrained>

%------------------------------- Reference --------------------------------
% J. Zhou, Y. Zhang, and P.N. Suganthan, “Constrained large-scale multiobjective optimization based on
% a competitive and cooperative swarm optimizer,” in Swarm and Evolutionary Computation, vol. 91, p.101735, 2024.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    methods
        function main(Algorithm,Problem)
            %% Generate random population
            [Rate,Acc] = Algorithm.ParameterSet(0.8,0.4);
            Population = Problem.Initialization();
            CV         = sum(max(0,Population.cons),2);
            CVmax      = max(CV);
            VAR0  = CVmax;
            if VAR0 == 0
                VAR0 = 1;
            end
            X=0;
            Competitive_Pop = UpdateP(Population,Problem.N,VAR0);
            [Cooperative_Pop,Cooperative_Pop_Fitness] = UpdateP(Population,Problem.N,inf);
            Population = ArchiveUpdate(Population,Problem.N);
            
            %% Optimization
            while Algorithm.NotTerminated(Population)
             %% Udate the epsilon value
                cp  = (-log(VAR0)-6)/log(1-0.5);
                VAR = VAR0*(1-X)^cp;
    
                
                Competitive_Pop_Fitness = CalFitness(Competitive_Pop.objs,Competitive_Pop.cons,VAR);
               
               [Off1Dec,Off1Vel] = CompetitiveOperator(Problem,Competitive_Pop,Competitive_Pop_Fitness);
                
                LearningPool = TournamentSelection(2,Problem.N,Cooperative_Pop_Fitness);
                Off2Dec    = CooperativeOperator(Problem,Cooperative_Pop(LearningPool).decs);
                
                  
                %% FDV
                if Problem.FE/Problem.maxFE <= Rate
                    Offspring1 = FDVOperator(Problem,Rate,Acc,Off1Dec,Off1Vel);
                    Offspring2 = FDVOperator(Problem,Rate,Acc,Off2Dec);
                else
                    Offspring1 = Problem.Evaluation(Off1Dec,Off1Vel);
                    Offspring2 = Problem.Evaluation(Off2Dec);
                end
                Offspring = [Offspring1,Offspring2];
                
                Population = ArchiveUpdate([Population,Offspring],Problem.N);
                Competitive_Pop = UpdateP([Competitive_Pop,Offspring],Problem.N,VAR);
                if Problem.FE/Problem.maxFE > Rate
                     [Cooperative_Pop,Cooperative_Pop_Fitness] = UpdateP([Offspring,Cooperative_Pop],Problem.N,VAR/Problem.N);
                else
                     [Cooperative_Pop,Cooperative_Pop_Fitness] = UpdateP([Offspring,Cooperative_Pop],Problem.N,inf);
                end

                 X = X + 1/(Problem.maxFE/Problem.N);
            end
        end
    end
end
