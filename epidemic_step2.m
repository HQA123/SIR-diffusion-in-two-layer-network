function new_states=epidemic_step2(old_states,graph,k,gamma)
    %���ɸ�Ⱦ������������ͬ�ĸ��崦�ڲ�ͬ��״̬�����ĸ�Ⱦ������ͬ,���崦�ڸ�Ⱦ̬���ʹ��ڸ�Ⱦ����������û�и�Ⱦ����
    infectiousness=(old_states+abs(old_states))/2;
    new_states=zeros(length(old_states),1);
    %��Ⱦ���������;�����ˣ���������Ч��Ⱦ�ʣ��γɸø�����һ�׶α���Ⱦ�ĸ���
    prob=(graph'*infectiousness)*k; %��Ⱦ��
    for individual=1:length(old_states)
         %���ڸ�Ⱦ̬�ĸ��壬�����һ�׶εĸ�Ⱦ̬�����߱������̬
        if (old_states(individual)>0)
            if (rand<infectiousness(individual)*gamma)
               new_states(individual)=-1; 
            else
               new_states(individual)=1;
            end
        else
             %ԭ��״̬���׸�Ⱦ̬���׸�Ⱦ̬��һ�����ʱ�Ϊ�Ѹ�Ⱦ̬
            if (old_states(individual)==0)
                if (rand<prob(individual))
                    new_states(individual)=1;
                else
                    new_states(individual)=0;

                end
%                 for i=1:length(old_states)
%                     if (infectiousness(i)==1)&(rand<graph(i,individual)*k)%����ط���΢�ַ��̵�k����һ��1/N��ת��
%                         new_states(individual)=1;
%                         break;
%                     else
%                         new_states(individual)=0
%                     end
%                 end
             %ԭ��������״̬�Ļ�������״̬
            else
                new_states(individual)=-1;
            end

        end
    end
end