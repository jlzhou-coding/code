classdef DPACMO < ALGORITHM
% <multi> <real/binary/permutation> <constrained>
% Coevolutionary constrained multi-objective optimization framework
% type --- 1 --- Type of operator (1. GA 2. DE)

%------------------------------- Reference --------------------------------
% J. Zhou, Y. Zhang and P.N. Suganthan, Dual population approximate constrained Pareto front for constrained 
% multiobjective optimization, Information Sciences, 119591, doi: https://doi.org/10.1016/j.ins.2023.119591.
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
            alpha=0.3;
            %% Generate random population
            Population1 = Problem.Initialization();
            Population2 = Problem.Initialization();
            
            cons = [Population1.cons; Population2.cons];
            epsilon0  = max(sum(max(0,cons),2));
            if epsilon0 == 0
                epsilon0 = 1;
            end
            r=0;
            cp=(-log(epsilon0)-6)/log(1-0.5);
            
            Fitness1    = CalFitness(Population1.objs,Population1.cons,0);
            Fitness2    = CalFitness(Population2.objs,Population2.cons,1e30);
            A     = ArchiveUpdate([Population1,Population2],Problem.N);
            %% Optimization
            while Algorithm.NotTerminated(Population1)
%                      Draw(Population2.objs,'sk','Markeredgecolor',[.5 .5 .5],'Markerfacecolor',[.8 .5 .5]);
                if Problem.FE/Problem.maxFE<=alpha
                     epsilon=1e30;
                else
                     epsilon = epsilon0.*(1-r).^cp;
                end
                if type == 1
                    MatingPool1 = TournamentSelection(2,Problem.N,Fitness1);
                    MatingPool2 = TournamentSelection(2,Problem.N,Fitness2);
                    Offspring1  = OperatorGAhalf(Population1(MatingPool1));
                    Offspring2  = OperatorGAhalf(Population2(MatingPool2));
                elseif type == 2
                    MatingPool1 = TournamentSelection(2,2*Problem.N,Fitness1);
                    MatingPool2 = TournamentSelection(2,2*Problem.N,Fitness2);
                    Offspring1  = OperatorDE(Population1,Population1(MatingPool1(1:end/2)),Population1(MatingPool1(end/2+1:end)));
                    Offspring2  = OperatorDE(Population2,Population2(MatingPool2(1:end/2)),Population2(MatingPool2(end/2+1:end)));

                end
                
                Q1=[Population1];
                Q2=[Population2];

                CV1 = sum(max(0,Q1.cons),2);
                Q1 = Q1(CV1<=epsilon);%只取可行解
                transfer_Population1= Q1;

                CV2 = sum(max(0,Q2.cons),2);
                Q2 = Q2(CV2<=0);%只取可行解
                transfer_Population2= Q2;

		 [Population1,Fitness1] = EnvironmentalSelection([Population1,Offspring1,Offspring2,transfer_Population2],Problem.N,0);
		 [Population2,Fitness2] = EnvironmentalSelection([Population2,Offspring2,Offspring1,transfer_Population1],Problem.N,epsilon);
                
               % Output the non-dominated and feasible solutions.
                A = ArchiveUpdate([A,Population1,Population2],Problem.N);
                if Problem.FE >= Problem.maxFE
                    Population1 = A;
                end
               r=Problem.FE/Problem.maxFE;
            end
        end
    end
end
