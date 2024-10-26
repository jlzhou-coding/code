classdef TNK < PROBLEM
% <problem> <Practical Problems>

%--------------------------------------------------------------------------
% M. Tanaka, H. Watanabe, Y. Furukawa, and T. Tanino, ¡°GA-based
% decision support system for multicriteria optimization,¡± in 1995 IEEE
% International Conference on Systems, Man and Cybernetics. Intelligent
% Systems for the 21st Century, vol. 2, Oct 1995, pp. 1556¨C1561 vol.2.
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
        function obj = TNK()
            obj.Global.M = 2;
            obj.Global.D = 2;
            if isempty(obj.Global.D)
                obj.Global.D = 2;
            end
            obj.Global.lower    = zeros(1,2)+10^(-20);
            obj.Global.upper    = pi*ones(1,2);
            obj.Global.encoding = 'real';
        end
        %% Calculate objective values
        function PopObj = CalObj(~,X)
           
            x = X';
            
            y(1,:) = x(1,:);
            y(2,:) = x(2,:);
            
            PopObj(:,1) = y(1,:)';
            PopObj(:,2) = y(2,:)';
        end
        %% Calculate constraint violations
        function PopCon = CalCon(~,X)
            x = X';
            c(1,:)= x(1,:) .^2 + x(2,:) .^2 - 1 - 0.1 * cos(16*atan(x(1,:) ./ x(2,:)));
            c(2,:)= 0.5 - (x(1,:) - 0.5) .^ 2 - (x(2,:) - 0.5) .^ 2 ;
            
            PopCon =-c';
        end
        %% Sample reference points on Pareto front
        function P = PF(~,~)
            f = textread('TNK.dat');
            P = {f(:,1:2)};
            P = P{1};
        end
    end
end
