%SIR ����ģ�����ؿ�����棬����ʱ�临�Ӷȣ�N*T��
function Epidemic_simulation_SIR()
%tt=zeros(101,31);
graph2=xlsread('C:\Users\HQA\Desktop\�����ָ�\����\�½��ļ���\���պͺ���granger2_0.01.xlsx'); %�����滻�ڽӾ���
graph1=xlsread('C:\Users\HQA\Desktop\�����ָ�\����\�½��ļ���\����ͺ���granger2_0.01.xlsx'); %�����
tic;
% sen_mat1 = zeros(11,11)
% sen_mat2 = zeros(11,11)
% for ibeta1 = 0:10
%     for ibeta2 = 0:10
        T=40;
        states_sum1=zeros(3,T);states_sum2=zeros(3,T);
        [N1,~]=size(graph1);[N2,~]=size(graph2);
%         beta1=0.0105;beta2=0.0440;
        beta1=0.0522; beta2=0.0019;
        mu1=0.0637; mu2=0.0811;
        %graph1(N1,:)=graph1(N1,:)*33;graph2(N2,:)=graph2 v(N2,:)*33; %����ý��ڵ�Ĵ�Ⱦ�������ӵ�N��
        for jj=1:100 %�ظ���ͼ��ƽ��
            old_states1=zeros(N1,1);old_states2=zeros(N2,1);
            % for ii=0:0.01:0.3
            %    N=5000;
            %     ��Ч������
            %    gph_count_max=3;%�м��Ⱦ״̬ͼ���
            %    global graph_count
            %    graph_count=1;
            %���߽ڵ�
            %     old_states([99,101,13,53,41,102],1)=-1;
            %     old_states([35,56,34,36,63,4,64,81,1],1)=-1;
            %     old_states([41,93,23,89,3,32,33,39,78],1)=-1;
            %     old_states(round(N*rand(1)),1)=8;
            
            %��ʼ����Դ
            %     nn=randperm(103,1)
            %     while ismember(nn,[99])
            %         nn=randperm(103,1)
            %     end
            %     old_states(nn,1)=1;
            %   ������A���ѡȡ�ڵ�Ϊ��ȾԴ
%             my_rand1 = randi([1 501]);my_rand2 = randi([1 501]); my_rand3 = randi([1 501]); 
            my_rand1 = randi([1 289]);my_rand2 = randi([1 289]); my_rand3 = randi([1 289]);
%             my_rand3 = randi([1 289]);my_rand4 = randi([1 289]);
            old_states1(my_rand1,1)=1; old_states1(my_rand2,1)=1; old_states1(my_rand3,1)=1; 
%             old_states2(my_rand3,1)=1; old_states2(my_rand4,1)=1;
            %     old_states1(128,1)=1;
            %     old_states1(87,1)=1;
            %     old_states1(181,1)=1;
            states1=zeros(N1,T);states2=zeros(N2,T);
            states_num1=zeros(3,T);states_num2=zeros(3,T);
            states_num1(:,1)=[N1-3 3 0]'; %��ʼ�У�ò��ûʲô�ã�
            states_num2(:,1)=[N2 0 0]'
            states1(:,1)=old_states1;
            states2(:,1)=old_states2;
            %     disease=[0 0 0.025 0.075 0.175 0.225 0.25 0.25];
            %     disease=[0.25 0.25 0.225 0.175 0.075 0.025 0 0];
            %     infectious_graph=zeros(N,N,3);
            for i=2:T
                %         states(:,i)=epidemic_step(states(:,i-1),graph,disease,k);
                states1(:,i)=epidemic_step2(states1(:,i-1),graph1,beta1,mu1);
%                 if states1(N1,i)==1
%                     disp(1)
%                 end
                states2(N2-49:N2,i-1)=states1(N1-49:N1,i);%��ͬ��ý��ڵ㻥��
                states2(:,i)=epidemic_step2(states2(:,i-1),graph2,beta2,mu2);
%                 if states2(N2,i)==1
%                     disp(2)
%                 end
                states1(N1-49:N1,i-1)=states2(N2-49:N2,i);%��ͬ��ý��ڵ㻥��,49Ϊý��Ľڵ���-1
                states_num1(:,i)=sir(states1(:,i)); %ͳ��SIR����
                states_num2(:,i)=sir(states2(:,i));
                %       if (i==3)||(i==6)||(i==10) %�׶θ�Ⱦͼ
                %            infectious_graph(:,:,graph_count-1)=graph_extract(states(:,i),graph,N);
                %       end
            end
            states_sum1=states_sum1+states_num1;
            states_sum2=states_sum2+states_num2;
        end %�ظ���ͼ��ƽ��
        states_sum1=states_sum1/100;states_sum2=states_sum2/100;
%         sen_mat1(ibeta1 + 1, ibeta2 + 1) = max(states_sum1(2, :)) / 551
%         sen_mat2(ibeta1 + 1, ibeta2 + 1) = max(states_sum2(2, :)) / 339
        %disp(max(states_sum1(2,:)));.
%     end
% end

figure(1);
h=plot(0:1:T-1,states_sum1(1,:),'g^-',0:1:T-1,states_sum1(2,:),'-sr',0:1:T-1,states_sum1(3,:),'-ob');
set(h,'LineWidth',2);
%     plot(T,sum((Y(dim*N+1:2*dim*N,:)-Y(1:dim*N,:))));
xlabel('time/day');ylabel('number of nodes');%title('��ʼ��Ⱦ��ҽҩ������ҵ');
set(gca,'FontSize',12);
set(get(gca,'XLabel'),'FontSize',18);%ͼ������Ϊ8 point��С5��
set(get(gca,'YLabel'),'FontSize',18);
set(get(gca,'title'),'FontSize',20);
legend('S','I','R');
axis([0 T 1 N1]);
%%%���ڶ���ͼ
figure(2);
h=plot(0:1:T-1,states_sum2(1,:),'g^-',0:1:T-1,states_sum2(2,:),'-sr',0:1:T-1,states_sum2(3,:),'-ob');
set(h,'LineWidth',2);
%     plot(T,sum((Y(dim*N+1:2*dim*N,:)-Y(1:dim*N,:))));
xlabel('time/day');ylabel('number of nodes');%title('��ʼ��Ⱦ��ҽҩ������ҵ');
set(gca,'FontSize',12);
set(get(gca,'XLabel'),'FontSize',18);%ͼ������Ϊ8 point��С5��
set(get(gca,'YLabel'),'FontSize',18);
set(get(gca,'title'),'FontSize',20);
legend('S','I','R');
axis([0 T 1 N2]);

%����м�׶θ�Ⱦͼ
%  for i=1:3
%      xlswrite(strcat('C:\Users\HQA\Desktop\����������ڷ��\����\�������\',num2str(i),'.xlsx'),infectious_graph(:,:,i));
%  end
toc;
%     tt(1,round(1+ii*100))=ii;
%     tt(jj+1,round(1+ii*100))=max(states_num(2,:));
% xlswrite('C:\Users\HQA\Desktop\�����ָ�\����\�½��ļ���\ph2��Ⱦϵ��������netA',sen_mat1);
% xlswrite('C:\Users\HQA\Desktop\�����ָ�\����\�½��ļ���\ph2��Ⱦϵ��������netB',sen_mat2);
end

function new_states=epidemic_step(old_states,graph,disease,k)
%���ɸ�Ⱦ������������ͬ�ĸ��崦�ڲ�ͬ��״̬�����ĸ�Ⱦ������ͬ
infectiousness=zeros(length(old_states),1);
global gamma;
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
        %             ���ո�Ⱦ������������
        if(old_states(individual)==length(disease)) %����Ⱦ�����־ͱ�Ϊ����
            new_states(individual)=-1;
        else
            new_states(individual)=old_states(individual)+1;
        end
        %             �����������ʷ���
        %             if (rand<=gamma)
        %                new_states(individual)=-1;
        %             end
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

function sir=sir(status); %ͳ��SIR������
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

%�м�ʱ���ĸ�Ⱦ�������˽ṹ
function A=graph_extract(state,graph,N)
A=zeros(N,N);
global graph_count
for j=1:N
    if state(j,1)>0
        A(j,:)=graph(j,:);
    end
end
graph_count=graph_count+1;
end