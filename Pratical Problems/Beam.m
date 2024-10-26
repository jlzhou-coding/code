classdef Beam < PROBLEM
% <problem> <Practical Problems>
% Beam

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
        function obj = Beam()
            if isempty(obj.Global.M)
                obj.Global.M =2;
            end
            if isempty(obj.Global.D)
                obj.Global.D = 4;
            end
            obj.Global.lower    = [0.125,0.1,0.1,0.125];
            obj.Global.upper    = [5.0,10,10,5.0];
            obj.Global.encoding = 'real';
        end
        %% Calculate objective values
        function PopObj = CalObj(~,PopDec)
%             M      = obj.Global.M;
            h=PopDec(:,1);
            l=PopDec(:,2);
            t=PopDec(:,3);
            b=PopDec(:,4);
            PopObj(:,1)=1.10471.*h.^2.*l+0.04811.*t.*b.*(14.0+t);
            PopObj(:,2)=2.1952./(b.*t.^3);
        end
        %% Calculate constraint violations
        function PopCon = CalCon(~,PopDec)
%             M      = obj.Global.M;
            h=PopDec(:,1);
            l=PopDec(:,2);
            t=PopDec(:,3);
            b=PopDec(:,4);
            tao_1=6000./(sqrt(2).*h.*l);
            tao_2 = (6000 .* (14.0 + 0.5 .* l) .* sqrt(0.25 .* (l.* l + (h + t) .* (h + t))))./ (2 .* (0.707 .* h .* l .* (l .* l ./ 12 + 0.25 .* (l .* l + (h + t) .* (h + t)))));
            theta=504000./(b.*t.^2);
            Pc=64764.022 .* (1 - 0.0282346 .* t) .* t .* b .* b .* b;
          
            PopCon(:,1) = -(13600 - sqrt(tao_1 .* tao_1 + tao_2 .* tao_2 + l .* tao_1 .* tao_2 ./ (sqrt(0.25 .* (l .* l + (h + t) .* (h + t))))));
		    PopCon(:,2) = -(30000 - theta);
		    PopCon(:,3) = -(b - h);
		    PopCon(:,4)= -(Pc - 6000);
        end
        %% A reference point for hypervolume calculation
        function P = PF(obj,N)
            CallStack = dbstack('-completenames');
            load(fullfile(fileparts(CallStack(1).file),'Beam_PF.mat'),'Beam_PF');
%             Fmax=max(Beam_PF);
%             Fmin=min(Beam_PF);
%             PF_normalized=(Beam_PF-Fmin)./(Fmax-Fmin);
             P = Beam_PF;
        end
    end
end