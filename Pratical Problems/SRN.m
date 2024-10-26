classdef SRN < PROBLEM
% <problem> <Practical Problems>

%--------------------------------------------------------------------------
% N. Srinivas and K. Deb, ¡°Multiobjective function optimization using
% nondominated sorting genetic algorithms,¡± IEEE Transactions on Evolutionary
% Computation, vol. 2, no. 3, pp. 1301¨C1308, 1994.
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
        function obj = SRN()
            obj.Global.M = 2;
            obj.Global.D = 2;
            if isempty(obj.Global.D)
                obj.Global.D = 2;
            end
            obj.Global.lower    = [-20 -20];
            obj.Global.upper    = [20 20];
            obj.Global.encoding = 'real';
        end
        %% Calculate objective values
        function PopObj = CalObj(~,X)
           
            x = X';
            
            y(1,:) = 2 + (x(1,:) - 2) .^ 2 + (x(2,:) - 1) .^ 2;
            y(2,:) = 9 * x(1,:) - (x(2,:) - 1) .^ 2;
            
            PopObj(:,1) = y(1,:)';
            PopObj(:,2) = y(2,:)';
        end
        %% Calculate constraint violations
        function PopCon = CalCon(~,X)
            x = X';
            c(1,:)= 225 - x(1,:) .^ 2 - x(2,:) .^ 2;
            c(2,:)= 3 * x(2,:) - x(1,:) - 10;
            
            PopCon =-c';
        end
        %% Sample reference points on Pareto front
        function P = PF(~,~)
            f = textread('SRN.dat');
            P = {f(:,1:2)};
            P = P{1};
        end
    end
end
