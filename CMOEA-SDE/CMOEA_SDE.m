function CMOEA_SDE(Global)
% <algorithm> <A>
% type --- 1 --- Type of operator (1. GA 2. DE)

%------------------------------- Reference --------------------------------
% J. Zhou, Y. Zhang, J. Zheng and M. Li, "Domination-Based Selection and Shift-Based Density Estimation for Constrained Multiobjective Optimization," 
% in IEEE Transactions on Evolutionary Computation, vol. 27, no. 4, pp. 993-1004, Aug. 2023, doi: 10.1109/TEVC.2022.3190401.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------
    
     %% Parameter setting
    type = Global.ParameterSet(1);
    
    %% Generate random population
    Population = Global.Initialization();
    CV=sum(max(Population.cons,0),2);
    r=Global.evaluated/Global.evaluation;
    Fitness = CalFitness(Population.objs,CV,r);
    arch = ArchiveUpdate(Population,Global.N);
   
    while Global.NotTermination(Population)
          r=Global.evaluated/Global.evaluation;
          if type == 1
            MatingPool = TournamentSelection(2,Global.N,Fitness);  
            Offspring  = GA(Population(MatingPool));
          elseif type == 2
            Mat1 = TournamentSelection(2,Global.N,Fitness);
            Mat2 = TournamentSelection(2,Global.N,Fitness);
            Offspring = DE(Population,Population(Mat1),Population(Mat2));
          end
        [Population,Fitness] = EnvironmentalSelection([Population,Offspring],Global.N,r);
		
		 % Output the non-dominated and feasible solutions.
         arch = ArchiveUpdate([arch,Population],Global.N);
         if Global.evaluated >= Global.evaluation
             Population = arch;
         end

    end
end