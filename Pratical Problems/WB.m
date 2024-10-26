classdef WB < PROBLEM
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
        function obj = WB()
            obj.Global.M = 2;
            obj.Global.D = 4;
            if isempty(obj.Global.D)
                obj.Global.D = 4;
            end
            obj.Global.lower    = [0.125 0.125 0.1 0.1];
            obj.Global.upper    = [5 5 10 10];
            obj.Global.encoding = 'real';
        end
        %% Calculate objective values
        function PopObj = CalObj(~,X)
           
            x = X';
            
            y(1,:) = 1.10471*x(1,:).^2.*x(3,:)+0.04811*x(4,:).*x(2,:).*(14+x(4,:));
            y(2,:) = 2.1952./(x(4,:).^3.*x(2,:));
            
            PopObj(:,1) = y(1,:)';
            PopObj(:,2) = y(2,:)';
        end
        %% Calculate constraint violations
        function PopCon = CalCon(~,X)
            x = X';
            Tao1 = 6000./(sqrt(2).*x(1,:).*x(3,:));
            Tao2 = 6000*(14+0.5*x(3,:)).*sqrt(0.25*(x(3,:).^2+(x(1,:)+x(4,:)).^2))./( 2*(0.707*x(1,:).*x(3,:).*( x(3,:).^2./12 + 0.25*(0.25*(x(3,:).^2+(x(1,:)+x(4,:)).^2))) ) );
            Tao = sqrt(Tao1.^2+Tao2.^2+ x(3,:).*Tao1.*Tao2./(sqrt(0.25*(x(3,:).^2+(x(1,:)+x(4,:)).^2))));
            Sigma  = 504000./(x(4,:).^2.*x(2,:));
            Pc = 64764.022*(1- 0.0282346*x(4,:)).*x(4,:).*x(2,:).^3;
            
            c(1,:) = 13600 - Tao;
            c(2,:) = 30000 - Sigma;
            c(3,:) = x(2,:) - x(1,:);
            c(4,:) = Pc - 6000;
            
            PopCon =-c';
        end
        %% Sample reference points on Pareto front
        function P = PF(~,~)
            f = textread('WB.dat');
            P = {f(:,1:2)};
            P = P{1};
        end
    end
end
