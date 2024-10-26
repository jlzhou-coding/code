classdef Brake < PROBLEM
% <problem> <Practical Problems>
% Brake

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
        function obj = Brake()
            if isempty(obj.Global.M)
                obj.Global.M =2;
            end
            if isempty(obj.Global.D)
                obj.Global.D = 4;
            end
            obj.Global.lower    = [55,75,1000,2];
            obj.Global.upper    = [80,110,3000,20];
            obj.Global.encoding = 'real';
        end
        %% Calculate objective values
        function PopObj = CalObj(~,PopDec)
%             M      = obj.Global.M;
            x1=PopDec(:,1);
            x2=PopDec(:,2);
            x3=PopDec(:,3);
            x4=PopDec(:,4);
            PopObj(:,1) = 4.9 .* 1e-5 .* (x2 .* x2 - x1 .* x1).*(x4 - 1);
		    PopObj(:,2) = 9.82 .* 1e6 .* (x2 .* x2 - x1 .* x1) ./ (x3 .* x4 .* (x2.^3 - x1.^3));
        end
        %% Calculate constraint violations
        function PopCon = CalCon(~,PopDec)
%             M      = obj.Global.M;
            x1=PopDec(:,1);
            x2=PopDec(:,2);
            x3=PopDec(:,3);
            x4=PopDec(:,4);
            PopCon(:,1) = -((x2 - x1) - 20);
            PopCon(:,2) = -(30 - 2.5 .* (x4 + 1));
            PopCon(:,3) = -(0.4 - x3 ./ (3.14 .* (x2 .* x2 - x1 .* x1)));
            PopCon(:,4) = -(1 - (2.22 .* 1e-3 .* x3 .* (x2.^3 - x1.^3)) ./(x2 .* x2 - x1 .* x1).^2);
            PopCon(:,5) = -((2.66 .* 1e-2 .* x3 .* x4 .* (x2.^3 - x1.^3)) ./ (x2 .* x2 - x1 .* x1) - 900);
        end
        %% A reference point for hypervolume calculation
        function P = PF(obj,N)
            CallStack = dbstack('-completenames');
            load(fullfile(fileparts(CallStack(1).file),'Brake_PF.mat'),'Brake_PF');
%             Fmax=max(Brake_PF);
%             Fmin=min(Brake_PF);
%             PF_normalized=(Brake_PF-Fmin)./(Fmax-Fmin);
             P = Brake_PF;
        end
    end
end