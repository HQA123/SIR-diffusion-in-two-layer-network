%SIR 传播模型蒙特卡洛仿真，运行时间复杂度：N*T。
function Epidemic_simulation_SIR();
    tic;
    N=5000;
    %有效传播率
    k=0.5;
    graph=BA_net(N);
    old_states=zeros(N,1);
    %old_states(round(N*rand(1)),1)=8;
    
    %初始化病源
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
    set(get(gca,'XLabel'),'FontSize',18);%图上文字为8 point或小5号
    set(get(gca,'YLabel'),'FontSize',18);
    axis([0 T 1 N]);
    hold on;
    toc;
   
end
function new_states=epidemic_step(old_states,graph,disease,k);
    %生成感染能力向量，因不同的个体处于不同的状态，它的感染能力不同
    infectiousness=zeros(length(old_states),1);
    for individual=1:length(old_states)
         %个体处于感染态，就存在感染能力，否者没有感染能力
        if (old_states(individual)>0)
            infectiousness(individual)=disease(old_states(individual));
        end  
    end
    %感染能力向量和矩阵相乘，并乘以有效感染率，形成该个体下一阶段被感染的概率
    prob=(graph*infectiousness)*k;
    for individual=1:length(old_states)
         %处于感染态的个体，变成下一阶段的感染态，或者变成免疫态
        if (old_states(individual)>0)
            if(old_states(individual)==length(disease))
                new_states(individual)=-1;
            else
                new_states(individual)=old_states(individual)+1;
            end
        else
             %原来状态是易感染态，易感染态以一定概率变为已感染态
            if (old_states(individual)==0)
                if (rand<prob(individual))
                    new_states(individual)=1;
                else
                    new_states(individual)=0;

                end
             %原来是免疫状态的还是免疫状态
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

%生成随即图
function A=suijitu(num)
    % N=input('网络图中节点的总数目N：');
    % alph=input('网络图中边的平均连接度alph:  ');
    % beta=input('表征边的平均长度的参数beta:  ');
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
%% 该函数求解两节点连接边的概率
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
 

%生成机制4：将概率归一化，利用赌轮法选择连接的边，直至生成边数为N*(N-1)/2*alph。
function A=bian_lianjie4(p,N,alph)
    A=zeros(N,N);
    p1=reshape(p,1,N*N)./sum(sum(p));   
    pp=cumsum(p1);%求累计概率
    k=0;
    while  k<N*(N-1)/2*alph           %利用赌轮法选择一条边相连
         random_data=rand(1,1);
         aa=find(pp>=random_data);jj=aa(1); % 节点jj即为用赌轮法选择的节点
         j=ceil(jj/N);i=jj-(j-1)*N;             %把单下标索引变为双下标索引，或者用函数ind2sub(siz,IND)
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

%生成无标度网络
function A=BA_net(num)
    %%% 从已有的m0个节点的网络开始，采用增长机制与优先连接的机制生成BA无标度网络
    %% A ——————返回生成网络的邻接矩阵
    % m0=input('未增长前的网络节点个数m0:  ');
    % m=input(' 每次引入的新节点时新生成的边数m： ');
    % N=input('增长后的网络规模N： ');
    % disp('初始网络时m0个节点的连接情况：1表示都是孤立；2表示构成完全图；3表示随机连接一些边');
    % pp=input('初始网络情况1，2或3： ');

    m0=15;
    m=5;
    N=num;

    pp=2;

    if m>m0
        disp('输入参数m不合法');
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
            disp('输入参数pp不合法');
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
        pp=cumsum(p);          %求累计概率
        for i=1:m              %利用赌轮法从已有的节点中随机选择m个节点与新加入的节点相连
            random_data=rand(1,1);
            aa=find(pp>=random_data);jj=aa(1); % 节点jj即为用赌轮法选择的节点
            A(k,jj)=-1;A(jj,k)=-1;
        end
    end
    n=size(A);
    for i=1:n
        A(i,i)=-sum(A(i,:));
    end

end


%生成小世界网络
function A=SW_net(num)
%%% 从有N个节点，每个节点有2K个邻居节点的最近邻耦合网络图通过随机化重连生成WS小世界网路
%% A ——————返回生成网络的邻接矩阵
% disp('该程序生成WS小世界网路：');
%请输入最近邻耦合网络中节点的总数N
N=num;
% K=input('请输入最近邻耦合网络中每个节点的邻居节点的个数的一半K：');
K=3;
% p=input('请输入随机化重连的概率p:');
p=0.2;
if K>floor(N/2)
    disp('输入的K值不合法')
    return;
end
% angle=0:2*pi/N:2*pi-2*pi/N;  %%生成最近邻耦合网络的各节点坐标
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
      A(i,jj)=-1; A(jj,i)=-1;     %%生成最近邻耦合网络的邻接矩阵
    end
end

for i=1:N
    for j=i+1:i+K
        jj=j;
        if j>N
            jj=mod(j,N);
        end
        p1=rand(1,1);
        if p1<p              %% 生成的随机数小于p，则边进行随机化重连,否则，边不进行重连
            A(i,jj)=0;A(jj,i)=0;  %重连策略：先断开原来的边，再在未连的边中随机选择另一个节点，与原节点连接。
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