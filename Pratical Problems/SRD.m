classdef SRD < PROBLEM
% <problem> <Practical Problems>

%--------------------------------------------------------------------------
% C. A. Coello Coello and G. T. Pulido, ¡°Multiobjective structural optimization
% using a microgenetic algorithm,¡± Structural and Multidisciplinary
% Optimization, vol. 30, no. 5, pp. 388¨C403, 2005.
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
        function obj = SRD()
            obj.Global.M = 2;
            obj.Global.D = 7;
            if isempty(obj.Global.D)
                obj.Global.D = 7;
            end
            obj.Global.lower    = [2.6 0.7 17 7.3 7.3 2.9 5.0];
            obj.Global.upper    = [3.6 0.8 28 8.3 8.3 3.9 5.5];
            obj.Global.encoding = 'real';
        end
        %% Calculate objective values
        function PopObj = CalObj(~,X)
           
            x = X';
            
            y(1,:) = 0.7854*x(1,:).*x(2,:).^2.*( 10*x(3,:).^2./3 + 14.933*x(3,:)- 43.0934 ) - 1.508*x(1,:).*( x(6,:).^2+x(7,:).^2 ) + 7.477*( x(6,:).^3+x(7,:).^3 )+0.7854*(x(4,:).*x(6,:).^2+x(5,:).*x(7,:).^2);
            y(2,:) = sqrt((745*x(4,:)./(x(2,:).*x(3,:))).^2 + 1.69*10.^7 )./(0.1*x(6,:).^3);
            
            PopObj(:,1) = y(1,:)';
            PopObj(:,2) = y(2,:)';
        end
        %% Calculate constraint violations
        function PopCon = CalCon(~,X)
            x = X';
            y(1,:) = 0.7854*x(1,:).*x(2,:).^2.*( 10*x(3,:).^2./3 + 14.933*x(3,:)- 43.0934 ) - 1.508*x(1,:).*( x(6,:).^2+x(7,:).^2 ) + 7.477*( x(6,:).^3+x(7,:).^3 )+0.7854*(x(4,:).*x(6,:).^2+x(5,:).*x(7,:).^2);
            y(2,:) = sqrt((745*x(4,:)./(x(2,:).*x(3,:))).^2 + 1.69*10.^7 )./(0.1*x(6,:).^3);
            c(1,:)= 1./(x(1,:).*x(2,:).^2.*x(3,:)) - 1./27;
            c(2,:)= 1./(x(1,:).*x(2,:).^2.*x(3,:).^2) - 1./397.5;
            c(3,:) = x(4,:).^3./(x(2,:).*x(3,:).*x(6,:).^4) - 1./1.93;
            c(4,:) = x(5,:).^3./(x(2,:).*x(3,:).*x(7,:).^4) - 1./1.93;
            c(5,:) = x(2,:).*x(3,:) - 40;
            c(6,:) = x(1,:)./x(2,:) - 12;
            c(7,:) = 5 - x(1,:)./x(2,:);
            c(8,:) = 1.9 - x(4,:) + 1.5*x(6,:);
            c(9,:) = 1.9 - x(5,:) + 1.1*x(7,:);
            c(10,:) = y(2,:) - 1300;
            c(11,:) = sqrt((745*x(5,:)./(x(2,:).*x(3,:))).^2 + 1.575*10.^8)./(0.1*x(7,:).^3)-1100;
            
            PopCon =c';
        end
        %% Sample reference points on Pareto front
        function P = PF(~,~)
            f = textread('SRD.dat');
            P = {f(:,1:2)};
            P = P{1};
        end
    end
end
