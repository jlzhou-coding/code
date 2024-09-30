function [Population,Fitness] = UpdateP(Population,N,epsilon)
% Selecte epsilon feasible solutions
            %% Calculate the fitness of each solution
             Fitness = CalFitness(Population.objs,Population.cons,epsilon);
            %% Non-dominated sorting
            CV=Population.cons;
            CV(CV<=epsilon)=0;
            [FrontNo,MaxFNo] = NDSort(Population.objs,CV,N);%Population.cons,
            Next = FrontNo <= MaxFNo;
            Population = Population(Next);
              Fitness=Fitness(Next);
            PCObj = Population.objs;
            nND   = length(Population);
            %% Population maintenance
            if length(Population) > N
                    % Normalization
                fmax  = max(PCObj,[],1);
                fmin  = min(PCObj,[],1);
                PCObj = (PCObj-repmat(fmin,nND,1))./repmat(fmax-fmin,nND,1);
                % Determine the radius of the niche
                d  = pdist2(PCObj,PCObj);
                d(logical(eye(length(d)))) = inf;
                sd = sort(d,2);
                r  = mean(sd(:,min(3,size(sd,2))));
                R  = min(d./r,1);
                % Delete solution one by one
                while length(Population) > N
                    [~,worst]  = max(1-prod(R,2));
                    Population(worst)  = [];
                    Fitness(worst)  = [];
                    R(worst,:) = [];
                    R(:,worst) = [];
                end
            end

end