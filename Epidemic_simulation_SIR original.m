%SIR ����ģ�����ؿ�����棬����ʱ�临�Ӷȣ�N*T��
function Epidemic_simulation_SIR();
    tic;
    N=5000;
    %��Ч������
    k=0.5;
    graph=BA_net(N);
    old_states=zeros(N,1);
    %old_states(round(N*rand(1)),1)=8;
    
    %��ʼ����Դ
    old_states(1:1,1)=4;
    T=200;
    states=zeros(N,T);
    states_num=zeros(3,T);
    states_num(:,1)=[N-1 1 0]';
    states(:,1)=old_states;
    disease=[0 0 0.025 0.075 0.175 0.225 0.25 0.25];
    for i=2:T
        states(:,i)=epidemic_step(states(:,i-1),graph,disease,k);
        states_num(:,i)=sir(states(:,i));
    end
    figure(1);
    h=plot(1:1:T,states_num(1,:),'g-',1:1:T,states_num(2,:),'r-',1:1:T,states_num(3,:),'b.');
    set(h,'LineWidth',2);
    % plot(T,sum((Y(dim*N+1:2*dim*N,:)-Y(1:dim*N,:))));
    xlabel('time');ylabel('number of individuals');
    set(gca,'FontSize',12); 
    set(get(gca,'XLabel'),'FontSize',18);%ͼ������Ϊ8 point��С5��
    set(get(gca,'YLabel'),'FontSize',18);
    axis([0 T 1 N]);
    hold on;
    toc;
   
end
function new_states=epidemic_step(old_states,graph,disease,k);
    %���ɸ�Ⱦ������������ͬ�ĸ��崦�ڲ�ͬ��״̬�����ĸ�Ⱦ������ͬ
    infectiousness=zeros(length(old_states),1);
    for individual=1:length(old_states)
         %���崦�ڸ�Ⱦ̬���ʹ��ڸ�Ⱦ����������û�и�Ⱦ����
        if (old_states(individual)>0)
            infectiousness(individual)=disease(old_states(individual));
        end  
    end
    %��Ⱦ���������;�����ˣ���������Ч��Ⱦ�ʣ��γɸø�����һ�׶α���Ⱦ�ĸ���
    prob=(graph*infectiousness)*k;
    for individual=1:length(old_states)
         %���ڸ�Ⱦ̬�ĸ��壬�����һ�׶εĸ�Ⱦ̬�����߱������̬
        if (old_states(individual)>0)
            if(old_states(individual)==length(disease))
                new_states(individual)=-1;
            else
                new_states(individual)=old_states(individual)+1;
            end
        else
             %ԭ��״̬���׸�Ⱦ̬���׸�Ⱦ̬��һ�����ʱ�Ϊ�Ѹ�Ⱦ̬
            if (old_states(individual)==0)
                if (rand<prob(individual))
                    new_states(individual)=1;
                else
                    new_states(individual)=0;

                end
             %ԭ��������״̬�Ļ�������״̬
            else
                new_states(individual)=-1;
            end

        end



    end
end

function history=epidemic(initstates,gragh,disease,k);

states=initstates;
count=sir(states);
history=count;
    while(count(2)>0)
        states=epidemic_step(states,graph,disease,k);
        count=sir(states);
        history(length(history(:,1))+1,:)=count;
    end

end

function sir=sir(status);
    sir=zeros(3,1);
    for individual=1:length(status);
       switch int32(status(individual))
           case -1
                sir(3)=sir(3)+1;
           case 0
               sir(1)=sir(1)+1;
               
           otherwise
               sir(2)=sir(2)+1;
       
       end
        
    end

end

%�����漴ͼ
function A=suijitu(num)
    % N=input('����ͼ�нڵ������ĿN��');
    % alph=input('����ͼ�бߵ�ƽ�����Ӷ�alph:  ');
    % beta=input('�����ߵ�ƽ�����ȵĲ���beta:  ');
    N=num;
    alph=0.02;
    beta=0.3;
    randData=rand(2,N)*1000;
    x=randData(1,:);
    y=randData(2,:);
    p=lianjiegailv(x,y,alph,beta,N);
    A=bian_lianjie4(p,N,alph);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% �ú���������ڵ����ӱߵĸ���
function p=lianjiegailv(x,y,alph,beta,N)
    d=zeros(N);
    for i=1:N
        for j=1:N
            d(i,j)=sqrt((x(i)-x(j))^2+((y(i)-y(j)))^2);
        end
    end
    L=max(max(d));
    for i=1:N
        for j=1:N
            p(i,j)=alph*exp(-d(i,j)/beta/L);
        end
            p(i,i)=0;
    end
end
 

%���ɻ���4�������ʹ�һ�������ö��ַ�ѡ�����ӵıߣ�ֱ�����ɱ���ΪN*(N-1)/2*alph��
function A=bian_lianjie4(p,N,alph)
    A=zeros(N,N);
    p1=reshape(p,1,N*N)./sum(sum(p));   
    pp=cumsum(p1);%���ۼƸ���
    k=0;
    while  k<N*(N-1)/2*alph           %���ö��ַ�ѡ��һ��������
         random_data=rand(1,1);
         aa=find(pp>=random_data);jj=aa(1); % �ڵ�jj��Ϊ�ö��ַ�ѡ��Ľڵ�
         j=ceil(jj/N);i=jj-(j-1)*N;             %�ѵ��±�������Ϊ˫�±������������ú���ind2sub(siz,IND)
        % [i,j=ind2sub(size(p),jj);
        if A(i,j)==0
            A(i,j)=-1;A(j,i)=-1;
            k=k+1;
        end   
    end
    n=size(A);
    for i=1:n
        A(i,i)=-sum(A(i,:));
    end
    
end

%�����ޱ������
function A=BA_net(num)
    %%% �����е�m0���ڵ�����翪ʼ�����������������������ӵĻ�������BA�ޱ������
    %% A ��������������������������ڽӾ���
    % m0=input('δ����ǰ������ڵ����m0:  ');
    % m=input(' ÿ��������½ڵ�ʱ�����ɵı���m�� ');
    % N=input('������������ģN�� ');
    % disp('��ʼ����ʱm0���ڵ�����������1��ʾ���ǹ�����2��ʾ������ȫͼ��3��ʾ�������һЩ��');
    % pp=input('��ʼ�������1��2��3�� ');

    m0=15;
    m=5;
    N=num;

    pp=2;

    if m>m0
        disp('�������m���Ϸ�');
        return;
    end
    x=100*rand(1,m0);
    y=100*rand(1,m0);

    switch  pp
        case 1
            A=zeros(m0);
        case 2
            A=-ones(m0);
            for i=1:m0
                A(i,i)=0;
            end
        case 3
            for i=1:m0
                for j=i+1:m0
                    p1=rand(1,1);
                    if p1>0.5 
                        A(i,j)=1;A(j,i)=0;
                    end
                end
            end
        otherwise
            disp('�������pp���Ϸ�');
            return;          
    end 
    for k=m0+1:N
        M=size(A,1);
        p=zeros(1,M);
        x0=100*rand(1,1);y0=100*rand(1,1);
        x(k)=x0;y(k)=y0;
        if isempty(find(A==-1))==0
            p(:)=1/M;
        else
             for i=1:M
                 p(i)=length(find(A(i,:)==-1))/length(find(A==-1));
             end
        end
        pp=cumsum(p);          %���ۼƸ���
        for i=1:m              %���ö��ַ������еĽڵ������ѡ��m���ڵ����¼���Ľڵ�����
            random_data=rand(1,1);
            aa=find(pp>=random_data);jj=aa(1); % �ڵ�jj��Ϊ�ö��ַ�ѡ��Ľڵ�
            A(k,jj)=-1;A(jj,k)=-1;
        end
    end
    n=size(A);
    for i=1:n
        A(i,i)=-sum(A(i,:));
    end

end


%����С��������
function A=SW_net(num)
%%% ����N���ڵ㣬ÿ���ڵ���2K���ھӽڵ��������������ͼͨ���������������WSС������·
%% A ��������������������������ڽӾ���
% disp('�ó�������WSС������·��');
%�������������������нڵ������N
N=num;
% K=input('��������������������ÿ���ڵ���ھӽڵ�ĸ�����һ��K��');
K=3;
% p=input('����������������ĸ���p:');
p=0.2;
if K>floor(N/2)
    disp('�����Kֵ���Ϸ�')
    return;
end
% angle=0:2*pi/N:2*pi-2*pi/N;  %%����������������ĸ��ڵ�����
% x=100*sin(angle);
% y=100*cos(angle);
% plot(x,y,'ro','MarkerEdgeColor','g','MarkerFaceColor','r','markersize',8);
% hold on; 

A=zeros(N);
for i=1:N
    for j=i+1:i+K
        jj=j;
        if j>N
            jj=mod(j,N);
        end
      A(i,jj)=-1; A(jj,i)=-1;     %%������������������ڽӾ���
    end
end

for i=1:N
    for j=i+1:i+K
        jj=j;
        if j>N
            jj=mod(j,N);
        end
        p1=rand(1,1);
        if p1<p              %% ���ɵ������С��p����߽������������,���򣬱߲���������
            A(i,jj)=0;A(jj,i)=0;  %�������ԣ��ȶϿ�ԭ���ıߣ�����δ���ı������ѡ����һ���ڵ㣬��ԭ�ڵ����ӡ�
            A(i,i)=inf; a=find(A(i,:)==0);
            rand_data=randint(1,1,[1,length(a)]);
            jjj=a(rand_data);
            A(i,jjj)=-1;A(jjj,i)=-1;
            A(i,i)=0;
        end
    end
end
A=-A;
%  for i=1:n
%         A(i,i)=-sum(A(i,:));
%  end

   
end