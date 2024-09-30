classdef SFEA < ALGORITHM
% <multi/many> <real/integer> <large/none> <constrained/none>
% type --- 1 --- Type of operator (1. GA 2. DE)

%------------------------------- Reference --------------------------------
% J. Zhou, Y. Zhang, X. Yang, Y. Fan, and P.N. Suganthan, “A staged fuzzy evolutionary algorithm for constrained large-scale multiobjective optimization,”
% in Applied Soft Computing, 2024.

%------------------------------- Copyright --------------------------------
% Copyright (c) 2022 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    methods
        function main(Algorithm,Problem)
            %% Parameter setting
            type = Algorithm.ParameterSet(1);

            %% Generate random population
            Population1 = Problem.Initialization();
            Population2 = Problem.Initialization();
            stage=0;
            Fitness1    = CalFitness(Population1.objs,Population1.cons,stage,true);
            Fitness2    = CalFitness(Population2.objs,Population2.cons,stage,false);
    
            %% Optimization
              arch = ArchiveUpdate([Population1,Population2],Problem.N);
            while Algorithm.NotTerminated(Population2)
                if type == 1
                    MatingPool1 = TournamentSelection(2,Problem.N,Fitness1);
                    MatingPool2 = TournamentSelection(2,Problem.N,Fitness2);
                    Off1Dec  = OperatorGAhalf(Problem,Population1(MatingPool1).decs);
                    Off1Dec = FDVOperator(Problem,Off1Dec,stage);
                    Offspring1 = Problem.Evaluation(Off1Dec);
                    Off2Dec  = OperatorGAhalf(Problem,Population2(MatingPool2).decs);
                    Off2Dec = FDVOperator(Problem,Off2Dec,stage);
                    Offspring2 = Problem.Evaluation(Off2Dec);
                elseif type == 2
                    MatingPool1 = TournamentSelection(2,2*Problem.N,Fitness1);
                    MatingPool2 = TournamentSelection(2,2*Problem.N,Fitness2);
                    Off1Dec  = OperatorDE(Problem,Population1.decs,Population1(MatingPool1(1:end/2)).decs,Population1(MatingPool1(end/2+1:end)).decs);
                    Off1Dec=FDVOperator(Problem,Off1Dec,stage);
                    Offspring1 = Problem.Evaluation(Off1Dec);
                    Off2Dec  = OperatorDE(Problem,Population2.decs,Population2(MatingPool2(1:end/2)).decs,Population2(MatingPool2(end/2+1:end)).decs);
                    Off2Dec=FDVOperator(Problem,Off2Dec,stage);
                    Offspring2 = Problem.Evaluation(Off2Dec);
                end
                [Population1,Fitness1] = EnvironmentalSelection([Population1,Offspring1,Offspring2],Problem.N,stage,true);
                [Population2,Fitness2] = EnvironmentalSelection([Population2,Offspring2,Offspring1],Problem.N,stage,false);
                % Output the non-dominated and feasible solutions.
                arch = ArchiveUpdate([arch,Population1,Population2],Problem.N);
                if Problem.FE >= Problem.maxFE
                    Population1 = arch;
                end
         
               r=Problem.FE/Problem.maxFE;
               sigmoid=16*1./(1+exp(-30*(r-0.7)));
               stage=round(sigmoid);
            end
        end
    end
end

