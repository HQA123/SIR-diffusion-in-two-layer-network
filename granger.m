
function [F,G,p]=granger(y,x,max_lag)
T = length(x);
BIC = zeros(max_lag,1);
RSSR = zeros(max_lag,1);
RSSR2 = zeros(max_lag,1);
i = 1;
while i <= max_lag
    ystar = x(i+1:T,:);
    xstar = [ones(T-i,1) zeros(T-i,i)];
    j = 1;
    while j <= i
        xstar(:,j+1) = x(i+1-j:T-j);
        j = j+1;
    end
    [~,~,r] = regress(ystar,xstar);
    BIC(i,:) = T*log(r'*r/T) + (i+1)*log(T);
    RSSR(i,:) = r'*r;
    RSSR2(i,:)=var(r); %计算残差的方差
    i = i+1;
end
%RSSR为不加滞后项x的受约束的残差平方和
x_lag = find(BIC==min(BIC));
r1=RSSR2(x_lag,:);
BIC = zeros(max_lag,1);
RSSUR = zeros(max_lag,1);
i = 1;
while i <= max_lag
    ystar = x(i+x_lag+1:T,:);
    xstar = [ones(T-(i+x_lag),1) zeros(T-(i+x_lag),x_lag+i)]; 
    j = 1;
    while j <= x_lag
        xstar(:,j+1) = x(i+x_lag+1-j:T-j,:);
        j = j+1;
    end
    %加入滞后项
    j = 1;
    while j <= i
        xstar(:,x_lag+j+1) = y(i+x_lag+1-j:T-j,:);
        j = j+1;
    end
    [~,bint,r] = regress(ystar,xstar);
    BIC(i,:) = T*log(r'*r/T) + (i+1)*log(T);
    RSSUR(i,:) = r'*r;
    RSSUR2(i,:)=var(r); %计算残差的方差
    i = i+1;
end
y_lag = find(BIC==min(BIC));
r2=RSSUR2(y_lag,:);
%感觉以下部分有些问题
F_num = ((RSSR(x_lag,:) - RSSUR(y_lag,:))/y_lag);
F_den = RSSUR(y_lag,:)/(T-(x_lag+y_lag+1));   
F = F_num/F_den;
%c_v = finv(1-alpha,y_lag,(T-(x_lag+y_lag+1)));
p = 1-fcdf(F,y_lag,(T-(x_lag+y_lag+1)));
G=log(r1/r2);
end