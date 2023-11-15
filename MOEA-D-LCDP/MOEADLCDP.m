classdef MOEADLCDP < ALGORITHM
    % <multi/many> <real/binary/permutation> <constrained/none>
    % MOEA/D based on LCDP, the simplest implementation.
    % type   ---   1 --- Type of operator (1. GA 2. DE)
    
    
    %------------------------------- Reference --------------------------------
    %J. Zhou, Y. Zhang, J. Wang and P. N. Suganthan, "Localized Constrained-Domination Principle for Constrained Multiobjective Optimization,"
	%in IEEE Transactions on Systems, Man, and Cybernetics: Systems, doi: 10.1109/TSMC.2023.3324797.

    %------------------------------- Copyright --------------------------------
    % Copyright (c) 2021 BIMK Group. You are free to use the PlatEMO for
    % research purposes. All publications which use this platform or any code
    % in the platform should acknowledge the use of "PlatEMO" and reference "Ye
    % Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
    % for evolutionary multi-objective optimization [educational forum], IEEE
    % Computational Intelligence Magazine, 2017, 12(4): 73-87".
    %--------------------------------------------------------------------------
    
    methods
        function main(Algorithm,Problem)
            %% Parameter setting
            type = Algorithm.ParameterSet(1);
            %% Generate the weight vectors
            [W,Problem.N] = UniformPoint(Problem.N,Problem.M);
            T = ceil(Problem.N/10);
            nr = ceil(Problem.N/100);
            
            %% Detect the neighbours of each solution
            B = pdist2(W,W);
            [~,B] = sort(B,2);
            B = B(:,1:T);
            
            %% Generate random population
            Population = Problem.Initialization();
            Population_ND = Population(NDSort(Population.objs,1)==1);
            z_max=max(Population_ND.objs,[],1);
            z_min=min(Population_ND.objs,[],1);
            Z = min(Population.objs,[],1);
            f_normalized=(Population.objs-repmat(z_min,size(Population,2),1))./(repmat(z_max,size(Population,2),1)-repmat(z_min,size(Population,2),1));
            f=sum(f_normalized,2);
            obj_mean(1)=mean(f);
            gen=1;
            %% Angle
            angle=acos(1-pdist2(W,W,'cosine'));
            temp_angle=angle;
            temp_angle(logical(eye(size(temp_angle))))=inf;
            theta_min=min(temp_angle');
            theta_min=theta_min';
            theta=theta_min.*0.5;
            
            arch    = ArchiveUpdate(Population,Problem.N);
            change=inf;
            %% Optimization
            while Algorithm.NotTerminated(Population)
                % For each solution
                for i = 1 : Problem.N
                    % Choose the parents
                    if rand < 0.9
                        P = B(i,randperm(size(B,2)));
                    else
                        P = randperm(Problem.N);
                    end
                    
                    % Generate an offspring
                    if type == 1
                        Offspring = OperatorGAhalf(Population(P(1:2)));
                    elseif type == 2
                        Offspring = OperatorDE(Population(i),Population(P(1)),Population(P(2)));
                    end
                    
                    % Update the ideal point
                    Z = min(Z,Offspring.obj);
                    
                    PopObj=Population(P).objs-repmat(Z,length(P),1);
              
                    Angle0   = acos(1 - pdist2(real(PopObj),W(P,:),'cosine'));
                    Angle = diag(Angle0);
                    NewAngle=acos(1-pdist2(real(Offspring.objs-Z),W(P,:),'cosine'));
                    NewAngle=NewAngle';
                
                    % Calculate the constraint violation of offspring and P
                    CVO = sum(max(0,Offspring.con));
                    CVP = sum(max(0,Population(P).cons),2);
                    
                    g_old = max(abs(Population(P).objs-repmat(Z,length(P),1))./W(P,:),[],2);
                    g_new = max(repmat(abs(Offspring.obj-Z),length(P),1)./W(P,:),[],2);
                    
                    eta=1e-6;
                    if  (change>=eta) 
                        case1=NewAngle<= theta(P) & Angle<= theta(P) & (g_old>=g_new & CVP==CVO| CVP>CVO);%both in niche
                        case2=NewAngle> theta(P) & Angle> theta(P) & (g_old>=g_new);%both not in niche
                        case3=NewAngle<= theta(P) & Angle> theta(P)& (g_old>=g_new);
                        case4=NewAngle> theta(P) & Angle<= theta(P) & g_old>=g_new & CVP>=CVO;
                        Population(P(find(case1 | case2 | case3 |case4,nr))) = Offspring;
                    else 
                        CVP=CVP.*(Angle+eta); CVO=CVO.*(NewAngle+eta);
                        Population(P(find(g_old>=g_new & CVP==CVO | CVP>CVO,nr))) = Offspring;
                    end
                end
                
                % Output the non-dominated and feasible solutions.
                arch = ArchiveUpdate([arch,Population],Problem.N);

                gen=gen+1;
                Population_ND = Population(NDSort(Population.objs,1)==1);
                z_max=max([Population_ND.objs;z_max],[],1);
                z_min=min([Population_ND.objs;z_min],[],1);
                f_normalized=(Population.objs-repmat(z_min,size(Population,2),1))./(repmat(z_max,size(Population,2),1)-repmat(z_min,size(Population,2),1));
                f=sum(f_normalized,2);
                obj_mean(gen)=mean(f);
                change=abs(obj_mean(gen)-obj_mean(gen-1));
                if Problem.FE >=Problem.maxFE
                    Population = arch;
                end
            end
        end
    end
end