classdef Two_Bar_Truss < PROBLEM
% <problem> <Practical Problems>
% Two Bar Truss

%------------------------------- Reference --------------------------------
% Multi Objective Design Optimization of two bar truss using NSGA II and TOPSIS
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
        function obj = Two_Bar_Truss()
            if isempty(obj.Global.M)
                obj.Global.M = 2;
            end
            if isempty(obj.Global.D)
                obj.Global.D = 2;
            end
            obj.Global.lower    = [0.10,2.00];%[0.10,2.25];
            obj.Global.upper    = [0.10,2.50];%[0.50,2.50];
            obj.Global.encoding = 'real';
        end
        %% Calculate objective values
        function PopObj = CalObj(~,PopDec)
%             M      = obj.Global.M;
            rou=0.283;
            h=100;
            P=10000;
            E=3.0E07;
            PopObj(:,1)=2*rou*h*PopDec(:,2).*sqrt(1+PopDec(:,1).^2);
            PopObj(:,2)=P*h*(1+PopDec(:,1).^2).^1.5 .* (1+PopDec(:,1).^4).^0.5./(2.*sqrt(2)*E.*PopDec(:,1).^2.*PopDec(:,2));
        end
        %% Calculate constraint violations
        function PopCon = CalCon(~,PopDec)
%             M      = obj.Global.M;
            P=10000;
            fai=2.0E04;
            PopCon(:,1)=P*(1+PopDec(:,1)).*(1+PopDec(:,1).^2).^0.5./(2.*sqrt(2).*PopDec(:,1).*PopDec(:,2))-fai;
            PopCon(:,2)=P*(1-PopDec(:,1)).*(1+PopDec(:,1).^2).^0.5./(2.*sqrt(2).*PopDec(:,1).*PopDec(:,2))-fai;
        end
         %% A reference point for hypervolume calculation
        function P = PF(obj,N)
              P = [200,1];
        end
    end
end