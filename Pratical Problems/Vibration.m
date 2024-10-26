classdef Vibration < PROBLEM
% <problem> <Practical Problems>
% Vibration Platform Design

%------------------------------- Reference --------------------------------
% T. Ray, K. Tai, and K. C. Seow, Multiobjective design optimization by an evolutionary algorithm, 
%Engineering Optimization, vol. 33, no. 4, pp. 399 424, 2001.
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
        function obj = Vibration()
            if isempty(obj.Global.M)
                obj.Global.M = 2;
            end
            if isempty(obj.Global.D)
                obj.Global.D = 5;
            end
            obj.Global.lower    = [3,0.35,0.01,0.01,0.01];
            obj.Global.upper    = [6,0.5,0.6,0.6,0.6];
            obj.Global.encoding = 'real';
        end
        %% Calculate objective values
        function PopObj = CalObj(~,PopDec)
%             M      = obj.Global.M;
            EI=(2.*PopDec(:,2)./3).*(70E09.*PopDec(:,3).^3+1.6E09.*(PopDec(:,4).^3-PopDec(:,3).^3)+200E09.*(PopDec(:,5).^3-PopDec(:,4).^3));
            mu=2.*PopDec(:,2).*(2770.*PopDec(:,3)+100.*(PopDec(:,4)-PopDec(:,3))+7780.*(PopDec(:,5)-PopDec(:,4)));
            PopObj(:,1)=-(pi./(2.*PopDec(:,1).^2)).*(EI./mu).^0.5;
            PopObj(:,2)=2.*PopDec(:,2).*(1500.*PopDec(:,3)+500.*(PopDec(:,4)-PopDec(:,3))+800.*(PopDec(:,5)-PopDec(:,4)));
        end
        %% Calculate constraint violations
        function PopCon = CalCon(~,PopDec)
%             M      = obj.Global.M;
            mu=2.*PopDec(:,2).*(2770.*PopDec(:,3)+100.*(PopDec(:,4)-PopDec(:,3))+7780.*(PopDec(:,5)-PopDec(:,4)));
            PopCon(:,1)=mu.*PopDec(:,1)-2800;
            PopCon(:,2)=0.00001-(PopDec(:,4)-PopDec(:,3));
            PopCon(:,3)=0.00001-(PopDec(:,5)-PopDec(:,4));
            PopCon(:,4)=PopDec(:,4)-PopDec(:,3)-0.01;
            PopCon(:,5)=PopDec(:,5)-PopDec(:,4)-0.01;
        end
        %% A reference point for hypervolume calculation
        function P = PF(obj,N)
              P = [0,300];
        end
    end
end