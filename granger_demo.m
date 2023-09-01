clc,clear
tic;
dat=xlsread('C:\Users\HQA\Desktop\外文B手稿\数据\新建文件夹\沪深和恒生return2.csv');%输入不带时间指标的矩阵
%dat=xlsread('C:\Users\HQA\Desktop\外文B手稿\数据\上证50和港.csv')
[nrow,ncol]=size(dat);
F=zeros(ncol,ncol);
sig=zeros(ncol,ncol);
for i=1:ncol
    for j=1:ncol
        if i==j
            continue;
        end
        [FF,G,p]=granger(dat(:,i),dat(:,j),3);%这里输入的回归阶数是最大回归阶数
        if p<=0.01  %置信系数
            sig(i,j)=1;
            F(i,j)=FF;
        end
    end
end
xlswrite('C:\Users\HQA\Desktop\外文B手稿\数据\新建文件夹\沪深和恒生granger2_0.01.xlsx',sig);
toc;
%xlswrite('C:\Users\HQA\Desktop\外文B手稿\数据\上证50和港.xlsx',F);