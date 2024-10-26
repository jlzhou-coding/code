classdef BNH < PROBLEM
% <problem> <Practical Problems>


%--------------------------------------------------------------------------
% T. T. Binh and U. Korn, ¡°MOBES: A multiobjective evolution strategy
% for constrained optimization problems,¡± in The Third International
% Conference on Genetic Algorithms (Mendel 97), vol. 25, 1997, p. 27.
%--------------------------------------------------------------------------

%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    methods
        %% Initialization
        function obj = BNH()
            obj.Global.M = 2;
            obj.Global.D = 2;
            if isempty(obj.Global.D)
                obj.Global.D = 2;
            end
            obj.Global.lower    = [0 0];
            obj.Global.upper    = [5 3];
            obj.Global.encoding = 'real';
        end
        %% Calculate objective values
        function PopObj = CalObj(~,X)
           
            x = X';
            y(1,:) = 4 * (x(1,:) .^ 2 + x(2,:) .^ 2) ;
            y(2,:) = (x(1,:) - 5) .^ 2 + (x(2,:) - 5) .^ 2;
            
            PopObj(:,1) = y(1,:)';
            PopObj(:,2) = y(2,:)';
        end
        %% Calculate constraint violations
        function PopCon = CalCon(~,X)
            x = X';
            c(1,:)= 25 - (x(1,:) - 5) .^ 2 - x(2,:) .^ 2;
            c(2,:)= (x(1,:) - 8) .^ 2 + (x(2,:) + 3) .^ 2 - 7.7;
            
            PopCon =-c';
        end
        %% Sample reference points on Pareto front
        function P = PF(~,~)
            f = textread('BNH.dat');
            P = {f(:,1:2)};
            P = P{1};
        end
    end
end
