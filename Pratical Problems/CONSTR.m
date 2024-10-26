classdef CONSTR < PROBLEM
% <problem> <Practical Problems>

%--------------------------------------------------------------------------
% P. D. Justesen, ¡°Multi-objective optimization using evolutionary algorithms,¡±
% University of Aarhus, Department of Computer Science,
% Denmark, 2009.
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
        function obj = CONSTR()
            obj.Global.M = 2;
            obj.Global.D = 2;
            if isempty(obj.Global.D)
                obj.Global.D = 2;
            end
            obj.Global.lower    = [0.1 0];
            obj.Global.upper    = [1  5];
            obj.Global.encoding = 'real';
        end
        %% Calculate objective values
        function PopObj = CalObj(~,X)
           
            x = X';
            y(1,:) = x(1,:);
            y(2,:) = (1 + x(2,:)) ./ x(1,:) ;
            
            PopObj(:,1) = y(1,:)';
            PopObj(:,2) = y(2,:)';
        end
        %% Calculate constraint violations
        function PopCon = CalCon(~,X)
            x = X';
            c(1,:)= x(2,:) + 9 * x(1,:) - 6;
            c(2,:)= -x(2,:) + 9 * x(1,:) - 1;
            
            PopCon =-c';
        end
        %% Sample reference points on Pareto front
        function P = PF(~,~)
            f = textread('CONSTR.dat');
            P = {f(:,1:2)};
            P = P{1};
        end
    end
end
