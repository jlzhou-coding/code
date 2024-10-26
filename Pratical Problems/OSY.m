classdef OSY < PROBLEM
% <problem> <Practical Problems>

%--------------------------------------------------------------------------
% A. Osyczka and S. Kundu, ¡°A new method to solve generalized
% multicriteria optimization problems using the simple genetic algorithm,¡±
% Structural Optimization, vol. 10, no. 2, pp. 94¨C99, Oct 1995.
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
        function obj = OSY()
            obj.Global.M = 2;
            obj.Global.D = 6;
            if isempty(obj.Global.D)
                obj.Global.D = 6;
            end
            obj.Global.lower    = [0 0 1 0 1 0];
            obj.Global.upper    = [10 10 5 6 5 10];
            obj.Global.encoding = 'real';
        end
        %% Calculate objective values
        function PopObj = CalObj(~,X)
           
            x = X';
            
            y(1,:) = -(25 * (x(1,:) - 2) .^ 2 + (x(2,:) -2) .^ 2 + (x(3,:) - 1) .^ 2 + (x(4,:) - 4) .^ 2 + (x(5,:) - 1) .^ 2);
            y(2,:) = sum(x .^ 2,1);
            
            PopObj(:,1) = y(1,:)';
            PopObj(:,2) = y(2,:)';
        end
        %% Calculate constraint violations
        function PopCon = CalCon(~,X)
            x = X';
            c(1,:) = x(1,:) + x(2,:) - 2;
            c(2,:) = 6 - x(1,:) - x(2,:);
            c(3,:) = 2 - x(2,:) + x(1,:);
            c(4,:) = 2 - x(1,:) + 3 * x(2,:);
            c(5,:) = 4 - (x(3,:) - 3) .^ 2 - x(4,:);
            c(6,:) = (x(5,:) - 3) .^ 2 + x(6,:) - 4;
            
            PopCon =-c';
        end
        %% Sample reference points on Pareto front
        function P = PF(~,~)
            f = textread('OSY.dat');
            P = {f(:,1:2)};
            P = P{1};
        end
    end
end
