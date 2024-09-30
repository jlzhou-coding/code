classdef MS_EA < ALGORITHM
% <multi> <real/integer/label/binary/permutation> <constrained>
% Constrained multiobjective evolutionary algorithm with multiple stages
% type   ---   2 --- Type of operator (1. GA 2. DE)
% lambda --- 0.5 --- Parameter for determining the current stage

%------------------------------- Reference --------------------------------
%周金龙,张英贵,肖杨,王娟. 不确定时间下多式联运多目标路径优化模型与算法[J/OL]. 交通运输系统工
%程与信息,2024. https://link.cnki.net/urlid/11.4520.u.20240913.1050.002
%------------------------------- Copyright --------------------------------
% Copyright (c) 2024 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    methods
        function main(Algorithm,Problem)
            %% Parameter setting
        	[type,lambda] = Algorithm.ParameterSet(2,0.5);

           %% Generate random population
            Population = Problem.Initialization();
            Fitness    = CalFitness(Population.objs);
            
            CV         = sum(max(0,Population.cons),2);
            CVmax      = max(CV);
            epsilon_0  = CVmax;
            epsilon    = epsilon_0;
            
            Tc    = 0.9 * ceil(Problem.maxFE/Problem.N);
            cp    = 2;
            alpha = 0.95;
            tao   = 0.05;

            %% Optimization
            while Algorithm.NotTerminated(Population)
                gen       = ceil(Problem.FE/Problem.N);
                CV        = sum(max(0,Population.cons),2);
                CV_max    = max(CV);
                CVmax     = max([CV_max,CVmax]);
                epsilon_0 = CVmax;
                rf = sum(CV <= 1e-6) / length(Population);
                epsilon = update_epsilon(tao,epsilon,epsilon_0,rf,alpha,gen,Tc,cp);
                
                if type == 1
                    MatingPool = TournamentSelection(2,Problem.N,Fitness);  
                    Offspring  = OperatorGA(Problem,Population(MatingPool));
                elseif type == 2
                    Mat1 = TournamentSelection(2,Problem.N,Fitness);
                    Mat2 = TournamentSelection(2,Problem.N,Fitness);
                    Offspring = OperatorDE(Problem,Population,Population(Mat1),Population(Mat2));
                end
                Q  = [Population,Offspring];
                  CV = sum(max(0,Q.cons),2);
                if mean(CV<=0) > lambda || Problem.FE >= lambda *Problem.maxFE %后期 All CV
                    [Population,Fitness] = UpdateP(Q,Problem.N,epsilon);% 多目标 epsilon
                else%前期
                      Fitness    = CalFitness(Q.objs);%收敛性、多样性
                      [Population,Fitness] = EnvironmentalSelection(Fitness,Q,Problem.N);
                end
                   A = ArchiveUpdate(Q,Problem.N);
                if Problem.FE >= Problem.maxFE
                    Population = A;
                end
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