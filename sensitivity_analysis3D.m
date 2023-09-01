clear all,clc;
graph4=xlsread('C:\Users\HQA\Desktop\外文手稿\数据\新建文件夹\敏感性分析\ph2传染系数敏感性100netB0.05_0.05_1.xls');
t=0:0.005:0.05;
[x,y]=meshgrid(t);%形成格点矩阵
surf(x,y,graph4)
surf(graph1);grid on
xlabel('x');
ylabel('y');
zlabel('f');
axis([0 0.05 0 0.05 0 0.8]);
figure(4)
surf(graph3)
title('f=cos[2*pi(2x-y)]; mesh')
% colormap winter
