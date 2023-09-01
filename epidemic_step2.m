function new_states=epidemic_step2(old_states,graph,k,gamma)
    %生成感染能力向量，因不同的个体处于不同的状态，它的感染能力不同,个体处于感染态，就存在感染能力，否者没有感染能力
    infectiousness=(old_states+abs(old_states))/2;
    new_states=zeros(length(old_states),1);
    %感染能力向量和矩阵相乘，并乘以有效感染率，形成该个体下一阶段被感染的概率
    prob=(graph'*infectiousness)*k; %感染力
    for individual=1:length(old_states)
         %处于感染态的个体，变成下一阶段的感染态，或者变成免疫态
        if (old_states(individual)>0)
            if (rand<infectiousness(individual)*gamma)
               new_states(individual)=-1; 
            else
               new_states(individual)=1;
            end
        else
             %原来状态是易感染态，易感染态以一定概率变为已感染态
            if (old_states(individual)==0)
                if (rand<prob(individual))
                    new_states(individual)=1;
                else
                    new_states(individual)=0;

                end
%                 for i=1:length(old_states)
%                     if (infectiousness(i)==1)&(rand<graph(i,individual)*k)%这个地方与微分方程的k存在一个1/N的转化
%                         new_states(individual)=1;
%                         break;
%                     else
%                         new_states(individual)=0
%                     end
%                 end
             %原来是免疫状态的还是免疫状态
            else
                new_states(individual)=-1;
            end

        end
    end
end