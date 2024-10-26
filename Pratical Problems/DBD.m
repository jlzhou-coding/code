classdef DBD < PROBLEM
% <problem> <Practical Problems>

%--------------------------------------------------------------------------
% T. Ray and K. M. Liew, ¡°A swarm metaphor for multiobjective design
% optimization,¡± Engineering Optimization, vol. 34, no. 2, pp. 141¨C153,
% 2002.
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
        function obj = DBD()
            obj.Global.M = 2;
            obj.Global.D = 4;
            if isempty(obj.Global.D)
                obj.Global.D = 4;
            end
            obj.Global.lower    = [55 75 1000 2];
            obj.Global.upper    = [80 110 3000 20];
            obj.Global.encoding = 'real';
        end
        %% Calculate objective values
        function PopObj = CalObj(~,X)
           
            x = X';
            y(1,:) = 4.9*10.^(-5).*(x(2,:).^2-x(1,:).^2).*(x(4,:) - 1);
            y(2,:) = 9.82*10.^6.*(x(2,:).^2-x(1,:).^2)./(x(3,:).*x(4,:).*(x(2,:).^3-x(1,:).^3));
            
            PopObj(:,1) = y(1,:)';
            PopObj(:,2) = y(2,:)';
        end
        %% Calculate constraint violations
        function PopCon = CalCon(~,X)
            x = X';
            c(1,:)= x(2,:) -x(1,:) - 20 ;
            c(2,:)= 30 - 2.5*(x(4,:) + 1);
            c(3,:) = 0.4-x(3,:)./(3.14*(x(2,:).^2-x(1,:).^2));
            c(4,:) = 1 - 2.22*10.^(-3).*x(3,:).*(x(2,:).^3-x(1,:).^3)./((x(2,:).^2-x(1,:).^2).^2);
            c(5,:) = 2.66*10.^(-2)*x(3,:).*x(4,:).*(x(2,:).^3-x(1,:).^3)./(x(2,:).^2-x(1,:).^2) - 900;
            PopCon =-c';
        end
        %% Sample reference points on Pareto front
        function P = PF(~,~)
            f = textread('DBD.dat');
            P = {f(:,1:2)};
            P = P{1};
        end
    end
end
