clear all,clc;
graph4=xlsread('C:\Users\HQA\Desktop\�����ָ�\����\�½��ļ���\�����Է���\ph2��Ⱦϵ��������100netB0.05_0.05_1.xls');
t=0:0.005:0.05;
[x,y]=meshgrid(t);%�γɸ�����
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
