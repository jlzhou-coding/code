classdef CarProblem < PROBLEM
% <problem> <Practical Problems>
% CarProblem

%------------------------------- Reference --------------------------------
% H. Jain and K. Deb, An evolutionary many-objective optimization algorithm
% using reference-point based non-dominated sorting approach, part II:
% Handling constraints and extending to an adaptive approach, IEEE
% Transactions on Evolutionary Computation, 2014, 18(4): 602-622.
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
        function obj = CarProblem()
            if isempty(obj.Global.M)
                obj.Global.M = 3;
            end
            if isempty(obj.Global.D)
                obj.Global.D = 7;
            end
            obj.Global.lower    = [0.5,0.45,0.5,0.5,0.875,0.4,0.4];
            obj.Global.upper    = [1.5,1.35,1.5,1.5,2.625,1.2,1.2];
            obj.Global.encoding = 'real';
        end
        %% Calculate objective values
        function PopObj = CalObj(~,PopDec)
%             M      = obj.Global.M;
            PopObj(:,1)=1.98+4.9*PopDec(:,1)+6.67*PopDec(:,2)+6.98*PopDec(:,3)+4.01*PopDec(:,4)+1.78*PopDec(:,5)+0.00001*PopDec(:,6)+2.73*PopDec(:,7);
            PopObj(:,2)=4.72-0.5*PopDec(:,4)-0.19*PopDec(:,2).*PopDec(:,3);
            PopObj(:,3)=0.5*(10.58-0.674*PopDec(:,1).*PopDec(:,2)-0.67275*PopDec(:,2)+16.45-0.489*PopDec(:,3).*PopDec(:,7)-0.843*PopDec(:,5).*PopDec(:,6));
        end
        %% Calculate constraint violations
        function PopCon = CalCon(~,PopDec)
%             M      = obj.Global.M;
            PopCon(:,1)=1.16-0.3717*PopDec(:,2).*PopDec(:,4)-0.0092928*PopDec(:,3)-1;
            PopCon(:,2)=0.261-0.0159*PopDec(:,1).*PopDec(:,2)-0.06486*PopDec(:,1)-0.019*PopDec(:,2).*PopDec(:,7)+0.0144*PopDec(:,3).*PopDec(:,5)+0.0154464*PopDec(:,6)-0.32;
            PopCon(:,3)=0.214+0.00817*PopDec(:,5)-0.045195*PopDec(:,1)-0.0135168*PopDec(:,1)+0.03099*PopDec(:,2).*PopDec(:,6)-0.018*PopDec(:,2).*PopDec(:,7)+0.007176*PopDec(:,3)+0.023232*PopDec(:,3)-0.00364*PopDec(:,5).*PopDec(:,6)-0.018*PopDec(:,2).^2-0.32;
            PopCon(:,4)=0.74-0.61*PopDec(:,2)-0.031296*PopDec(:,3)-0.031872*PopDec(:,7)+0.227*PopDec(:,2).^2-0.32;
            PopCon(:,5)=28.98+3.818*PopDec(:,3)-4.2*PopDec(:,1).*PopDec(:,2)+1.27296*PopDec(:,6)-2.68065*PopDec(:,7)-32;
            PopCon(:,6)=33.86+2.95*PopDec(:,3)-5.057*PopDec(:,1).*PopDec(:,2)-3.795*PopDec(:,2)-3.4431*PopDec(:,7)+1.45728-32;
            PopCon(:,7)=46.36-9.9*PopDec(:,2)-4.4505*PopDec(:,1)-32;
            PopCon(:,8)=4.72-0.5*PopDec(:,4)-0.19*PopDec(:,2).*PopDec(:,3)-4;
            PopCon(:,9)=10.58-0.674*PopDec(:,1).*PopDec(:,2)-0.67275*PopDec(:,2)-9.9;
            PopCon(:,10)=16.45-0.489*PopDec(:,3).*PopDec(:,7)-0.843*PopDec(:,5).*PopDec(:,6)-15.7;
        end
        %% Sample reference points on Pareto front
        function P = PF(obj,N)
            P=[50,5,15];
        end
    end
end