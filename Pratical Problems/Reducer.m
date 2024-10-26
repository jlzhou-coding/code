classdef Reducer < PROBLEM
% <problem> <Practical Problems>
% Reducer

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
        function obj = Reducer()
            if isempty(obj.Global.M)
                obj.Global.M =2;
            end
            if isempty(obj.Global.D)
                obj.Global.D = 7;
            end
            obj.Global.lower    = [2.6,0.7,17.0,7.3,7.3,2.9,5.0];
            obj.Global.upper    = [3.6,0.8,28.0,8.3,8.3,3.9,5.5];
            obj.Global.encoding = 'real';
        end
        %% Calculate objective values
        function PopObj = CalObj(~,PopDec)
%             M      = obj.Global.M;
            PopDec(:,3)=round(PopDec(:,3));
            sum1 = 0.7854 .* PopDec(:,1) .* PopDec(:,2) .* PopDec(:,2) .* (10.0 .* PopDec(3) .* PopDec(:,3) / 3 + 14.933 .* PopDec(:,3) - 43.0934);
            sum2 = 1.508 .* PopDec(:,1) .* (PopDec(:,6) .* PopDec(:,6) + PopDec(:,7) .* PopDec(:,7)) - 7.477 .* (PopDec(:,6).^3 + PopDec(:,7).^3);
            sum3 = 0.7854 .* (PopDec(:,4) .* PopDec(:,6) .* PopDec(:,6) + PopDec(:,5) .* PopDec(:,7) .* PopDec(:,7));
            PopObj(:,1)=sum1-sum2+sum3;
            PopObj(:,2)=sqrt((745 .* PopDec(:,4)./(PopDec(:,2).* PopDec(:,3))).^2 + 1.69.* 1e7)./ (0.1.* PopDec(:,6).^3);
        end
        %% Calculate constraint violations
        function PopCon = CalCon(~,PopDec)
%             M      = obj.Global.M;
            PopObj(:,2)=sqrt((745 .* PopDec(:,4)./(PopDec(:,2).* PopDec(:,3))).^2 + 1.69.* 1e7)./ (0.1.* PopDec(:,6).^3);
            PopCon(:,1) = -(1.0 ./ 27 - 1.0 ./ (PopDec(:,1) .* PopDec(:,2) .* PopDec(:,2) .*PopDec(:,3)));
            PopCon(:,2) = -(1.0 ./ 397.5 - 1.0 ./ (PopDec(:,1) .* PopDec(:,2) .* PopDec(:,2) .* PopDec(:,3) .* PopDec(:,3)));
            PopCon(:,3) = -(1 ./ 1.93 - PopDec(:,4).^3 ./ (PopDec(:,2) .* PopDec(:,3) .* (PopDec(:,6).^4)));
            PopCon(:,4) = -(1 ./ 1.93 - PopDec(:,5).^3 ./ (PopDec(:,2) .* PopDec(:,3) .* PopDec(:,7).^4));
            PopCon(:,5) = -(40 - PopDec(:,2) .* PopDec(:,3));
            PopCon(:,6) = -(12 - PopDec(:,1) ./ PopDec(:,2));
            PopCon(:,7) = -(PopDec(:,1) ./ PopDec(:,2) - 5);
            PopCon(:,8) = -(PopDec(:,4) - 1.5 .* PopDec(:,6) - 1.9);
            PopCon(:,9) = -(PopDec(:,5) - 1.1 .* PopDec(:,7) - 1.9);
            PopCon(:,10) = -(1300 - PopObj(:,2));
            PopCon(:,11) = -(1100 - sqrt((745 .* PopDec(:,5) ./ (PopDec(:,2) .* PopDec(:,3))).^2 + 1.575 .* 1e8) ./ (0.1 .* PopDec(:,7).^3));
        end
        %% A reference point for hypervolume calculation
        function P = PF(obj,N)
            CallStack = dbstack('-completenames');
            load(fullfile(fileparts(CallStack(1).file),'Reducer_PF.mat'),'Reducer_PF');
%             Fmax=max(Reducer_PF);
%             Fmin=min(Reducer_PF);
%             PF_normalized=(Reducer_PF-Fmin)./(Fmax-Fmin);
             P = Reducer_PF;
        end
    end
end