clc,clear
tic;
dat=xlsread('C:\Users\HQA\Desktop\����B�ָ�\����\�½��ļ���\����ͺ���return2.csv');%���벻��ʱ��ָ��ľ���
%dat=xlsread('C:\Users\HQA\Desktop\����B�ָ�\����\��֤50�͸�.csv')
[nrow,ncol]=size(dat);
F=zeros(ncol,ncol);
sig=zeros(ncol,ncol);
for i=1:ncol
    for j=1:ncol
        if i==j
            continue;
        end
        [FF,G,p]=granger(dat(:,i),dat(:,j),3);%��������Ļع���������ع����
        if p<=0.01  %����ϵ��
            sig(i,j)=1;
            F(i,j)=FF;
        end
    end
end
xlswrite('C:\Users\HQA\Desktop\����B�ָ�\����\�½��ļ���\����ͺ���granger2_0.01.xlsx',sig);
toc;
%xlswrite('C:\Users\HQA\Desktop\����B�ָ�\����\��֤50�͸�.xlsx',F);