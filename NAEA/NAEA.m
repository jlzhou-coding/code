function NAEA(Global)
% <algorithm> <N>
%------------------------------- Reference --------------------------------
% J. Zhou, J. Zou, S. Yang, J. Zheng, D. Gong, and T. Pei, ¡°Niche-based
%and angle-based selection strategies for many-objective evolutionary
%optimization,¡± Information Sciences, vol. 571, pp. 133¨C153, 2021.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    %% Generate random population
    Population = Global.Initialization();
     FrontNo    = NDSort(Population.objs,Population.cons,inf);
    %% Optimization
    while Global.NotTermination(Population)
       % MatingPool = randi(Global.N,1,Global.N);
	    MatingPool = TournamentSelection(2,Global.N,FrontNo,sum(max(0,Population.cons),2));
        Offspring  = GA(Population(MatingPool));  
        Population = [Population,Offspring];
        [Population,FrontNo] = EnvironmentalSelection(Population,Global.N);
    end
end