function MOEAD_DMEpsilon(Global)
% <algorithm> <C>
% An improved epsilon constraint-handling method in MOEA/D
% for CMOPs with large infeasible regions
%  Z. Fan, H. Li, C. Wei, W. Li, H. Huang, X. Cai, and Z. Cai, “An
% improved epsilon constraint handling method embedded in MOEA/D for
% constrained multi-objective optimization problems,” in Computational
% Intelligence (SSCI), 2016 IEEE Symposium Series on. IEEE, 2016, pp.
% 1–8.
% type --- 1 --- The type of aggregation function

%------------------------------- Reference --------------------------------
%  M. Asafuddoula, T. Ray, R. Sarker, and K. Alam, “An adaptive constraint
%handling approach embedded MOEA/D,” in Evolutionary Computation
%(CEC), 2012 IEEE Congress on. IEEE, 2012, pp. 1–8.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------
   
    %% Parameter setting
    type = Global.ParameterSet(2);
 
    %% Generate the weight vectors
    [W,Global.N] = UniformPoint(Global.N,Global.M);
        T=30;
        nr=2;
    %% Detect the neighbours of each solution
    B = pdist2(W,W);
    [~,B] = sort(B,2);
    B = B(:,1:T);
    
    %% Generate random population
    Population = Global.Initialization();
    Z = min(Population.objs,[],1);
    %% Evaluate the Population
    alpha=0.95;
    tao=0.1;
    CVI        = sum(max(0,Population.cons),2);
    phi_max=max(CVI);
    epsilon_0=(1+alpha)*max(CVI);
    epsilon_k=epsilon_0;
    arch = ArchiveUpdate(Population,Global.N);
    c=1;
 
    delta=0.8;%alpha%%%%%%%%0.8
    %% Optimization
    while Global.NotTermination(Population)
         CV        = sum(max(0,Population.cons),2);
         phi_max=max([phi_max;CV]);
        infindex = (CV ~= 0);%infeasible solutions
        rf = 1 - sum(infindex)/length(CV);%rate of feasible
 
       t=1*Global.N;%改变theta0,rf小，N小，theta大
        if Global.gen<delta*Global.maxgen
%         theta = (2/Global.N/2)*(1+Global.gen/Global.maxgen)^(log(Global.N)/log(1+delta));
          theta = (2/t/2)*(1+Global.gen/Global.maxgen)^(log(t)/log(1+delta));
        else
            theta = 1;
        end
        x(c)= theta;
        c=c+1;
   

              epsilon_k =  update_epsilon(tao,epsilon_k,epsilon_0,rf,alpha,phi_max,theta,Global.gen,0.8 * Global.maxgen);   



        % For each solution
        for i = 1 : Global.N      
           % Choose the parents
            if rand < 0.90
                P = B(i,randperm(size(B,2)));
            else
                P = randperm(Global.N);
            end

            % Generate an offspring
            if contains(class(Global.problem),'LIRCMOP') || contains(class(Global.problem),'DOC')
                Offspring = DE(Population(i),Population(P(1)),Population(P(2)));
            else
                Offspring = GAhalf(Population(P(1:2)));
            end

            % Update the ideal point
            Z = min(Z,Offspring.obj);
			
			% Calculate the constraint violation of offspring and P
            CVO = sum(max(0,Offspring.con));
            CVP = sum(max(0,Population(P).cons),2);
     
         
            % Update the neighbours
            switch type
                case 1
                    % PBI approach
                    normW   = sqrt(sum(W(P,:).^2,2));
                    normP   = sqrt(sum((Population(P).objs-repmat(Z,length(P),1)).^2,2));
                    normO   = sqrt(sum((Offspring.obj-Z).^2,2));
                    CosineP = sum((Population(P).objs-repmat(Z,length(P),1)).*W(P,:),2)./normW./normP;
                    CosineO = sum(repmat(Offspring.obj-Z,length(P),1).*W(P,:),2)./normW./normO;
                    g_old   = normP.*CosineP + 5*normP.*sqrt(1-CosineP.^2);
                    g_new   = normO.*CosineO + 5*normO.*sqrt(1-CosineO.^2);
                case 2
                    % Tchebycheff approach
                    g_old = max(abs(Population(P).objs-repmat(Z,length(P),1)).*W(P,:),[],2);
                    g_new = max(repmat(abs(Offspring.obj-Z),length(P),1).*W(P,:),[],2);
                case 3
                    % Tchebycheff approach with normalization
                    Zmax  = max(Population.objs,[],1);
                    g_old = max(abs(Population(P).objs-repmat(Z,length(P),1))./repmat(Zmax-Z,length(P),1).*W(P,:),[],2);
                    g_new = max(repmat(abs(Offspring.obj-Z)./(Zmax-Z),length(P),1).*W(P,:),[],2);
                case 4
                    % Modified Tchebycheff approach
                    g_old = max(abs(Population(P).objs-repmat(Z,length(P),1))./W(P,:),[],2);
                    g_new = max(repmat(abs(Offspring.obj-Z),length(P),1)./W(P,:),[],2);
            end
           CV_max=max(CVO,CVP);
           g_max=max(g_old,g_new);
       fit_old=(1-theta)*(g_old./(g_max+eps))+(theta)*(CVP./(CV_max+eps));
       fit_new=(1-theta)*(g_new./(g_max+eps))+(theta)*(CVO./(CV_max+eps));
        
       
             Population(P(find(((g_old>=g_new) & (CVP==CVO | (CVP<=epsilon_k & CVO<=epsilon_k))) | fit_old>fit_new,nr))) = Offspring;
          
        end
        % Output the non-dominated and feasible solutions.
        arch = ArchiveUpdate([arch,Population],Global.N);
        if Global.evaluated >= Global.evaluation
            Population = arch;
            OBJ=Population.objs;
            y=0;
        end
    end
end

function result = update_epsilon(tao,epsilon_k,epsilon_0,rf,alpha,phi_max,theta,gen,Tc)%theta in [0,1]
if rf < alpha
    result = (1 - tao) * epsilon_k;
else
    result =epsilon_0;
end

end