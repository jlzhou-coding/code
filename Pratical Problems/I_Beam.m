classdef I_Beam < PROBLEM
% <problem> <Practical Problems>
% I_Beam

%------------------------------- Reference --------------------------------
% MOEA/D with angle-based constrained dominance principle for constrained multi-objective optimization problems 
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
        function obj = I_Beam()
            if isempty(obj.Global.M)
                obj.Global.M = 2;
            end
            if isempty(obj.Global.D)
                obj.Global.D = 4;
            end
            obj.Global.lower    = [10,10,0.9,0.9];
            obj.Global.upper    = [80,50,5,5];
            obj.Global.encoding = 'real';
        end
        %% Calculate objective values
        function PopObj = CalObj(~,PopDec)
%             M      = obj.Global.M;
            E=20000;
            P=600;
            l=200;
            I=(PopDec(:,3).*(PopDec(:,1)-2.*PopDec(:,4)).^3+2.*PopDec(:,2).*PopDec(:,4).*(4.*PopDec(:,4).^2+3.*PopDec(:,1).*(PopDec(:,1)-2.*PopDec(:,4))))./12;
           
            PopObj(:,1)=2.*PopDec(:,2).*PopDec(:,4)+PopDec(:,3).*(PopDec(:,1)-2.*PopDec(:,4));
            PopObj(:,2)=P.*l.^3./(48.*E.*I);
        end
        %% Calculate constraint violations
        function PopCon = CalCon(~,PopDec)
%             M      = obj.Global.M;
            kg=16;
            Wy=(PopDec(:,3).*(PopDec(:,1)-2.*PopDec(:,4)).^3+2.*PopDec(:,2).*PopDec(:,4).*(4.*PopDec(:,4).^2+3.*PopDec(:,1).*(PopDec(:,1)-2.*PopDec(:,4))))./(6.*PopDec(:,1));
            Wz=((PopDec(:,1)-2.*PopDec(:,4)).*PopDec(:,3).^3+2.*PopDec(:,4).*PopDec(:,2).^3)./(6.*PopDec(:,2));
            My=30000;
            Mz=2500;
            PopCon(:,1)=-1*(kg-My./Wy-Mz./Wz);
        end
        %% A reference point for hypervolume calculation
        function P = PF(obj,N)
            CallStack = dbstack('-completenames');
            load(fullfile(fileparts(CallStack(1).file),'Ibeam_PF.mat'),'Ibeam_PF');
%             Fmax=max(Ibeam_PF);
%             Fmin=min(Ibeam_PF);
%             PF_normalized=(Ibeam_PF-Fmin)./(Fmax-Fmin);
            P = Ibeam_PF;
        end
    end
end