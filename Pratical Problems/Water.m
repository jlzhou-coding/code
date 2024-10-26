classdef Water < PROBLEM
% <problem> <Practical Problems>
% Water Problem

%------------------------------- Reference --------------------------------
%K. Musselman and J. Talavage, A tradeoff cut approach to mul- tiple objective
% optimization, Operations Research, vol. 28, no. 6, pp. 1424 1435, 1980.
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
        function obj = Water()
            if isempty(obj.Global.M)
                obj.Global.M =5;
            end
            if isempty(obj.Global.D)
                obj.Global.D = 3;
            end
            obj.Global.lower    = [0.01,0.01,0.01];
            obj.Global.upper    = [0.45,0.10,0.10];
            obj.Global.encoding = 'real';
        end
        %% Calculate objective values
        function PopObj = CalObj(~,PopDec)
%             M      = obj.Global.M;
            PopObj(:,1)=106780.37*(PopDec(:,2)+PopDec(:,3))+61704.67;
            PopObj(:,2)=3000.0*PopDec(:,1);
            PopObj(:,3)=305700*2289*PopDec(:,2)./(0.06*2289.0).^0.65;
            PopObj(:,4)=250*2289*exp(-1*39.75*PopDec(:,2)+9.9*PopDec(:,3)+2.74);
            PopObj(:,5)=25*(1.39./(PopDec(:,1).*PopDec(:,2))+4940*PopDec(:,3)-80);
        end
        %% Calculate constraint violations
        function PopCon = CalCon(~,PopDec)
%             M      = obj.Global.M;
            PopCon(:,1)=0.00139./(PopDec(:,1).*PopDec(:,2))+4.94*PopDec(:,3)-0.08-1;
            PopCon(:,2)=0.000306./(PopDec(:,1).*PopDec(:,2))+1.082*PopDec(:,3)-0.0986-1;
            PopCon(:,3)=12.307./(PopDec(:,1).*PopDec(:,2))+49408.24*PopDec(:,3)+4051.02-50000;
            PopCon(:,4)=2.098./(PopDec(:,1).*PopDec(:,2))+8046.33*PopDec(:,3)-696.71-16000;
            PopCon(:,5)=2.138./(PopDec(:,1).*PopDec(:,2))+7883.39*PopDec(:,3)-705.04-10000;
            PopCon(:,6)=0.417.*(PopDec(:,1).*PopDec(:,2))+1721.26*PopDec(:,3)-136.54-2000;
            PopCon(:,7)=0.164./(PopDec(:,1).*PopDec(:,2))+631.13*PopDec(:,3)-54.48-550;
        end
        %% A reference point for hypervolume calculation
        function P = PF(obj,N)
             CallStack = dbstack('-completenames');
             load(fullfile(fileparts(CallStack(1).file),'Water_PF.mat'),'Water_PF');
%              Fmax=max(Water_PF);
%              Fmin=min(Water_PF);
%              PF_normalized=(Water_PF-Fmin)./(Fmax-Fmin);
             P = Water_PF;
        end
    end
end