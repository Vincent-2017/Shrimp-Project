clear ;clc ;close all ;
% 二值化
img = imread('8.jpg');
[h,w,n] = size(img);
hsv = rgb2hsv(img);  % 转hsv
s = hsv(:,:,2);      % s分量
bw = ~im2bw(s,graythresh(s)); % 二值化
se = strel('disk',10);
bw2 = imclose(bw,se);
bw3 = bwareaopen(bw2,500); % 二值化的最终结果

% 求边缘
bbw3 = double(bw3);
bwedge = edge(bbw3);
indexedge = find(bwedge ~= 0);
[rowedge,coledge] = ind2sub(size(bwedge),indexedge);
Yedge = rowedge';
Xedge = coledge';
lenedge = length(Yedge);

% 求边缘图像的质心
Xsum = 0;
Ysum = 0;
for i = 1:lenedge
    Xsum = Xsum + Xedge(i);
    Ysum = Ysum + Yedge(i);
end
Xcenter = Xsum / lenedge ;
Ycenter = Ysum / lenedge ; % 质心

% matlab细化
bw4 = bwmorph(bw3,'thin',Inf); % 原始细化图
indexbw4 = find(bw4 ~= 0); % 非零点坐标索引
[rowbw4,colbw4] = ind2sub(size(bw4),indexbw4);
Ybw4 = rowbw4';
Xbw4 = colbw4';
lenbw4 = length(Xbw4);

% 显示
showimg = bwedge;
showimg(bw4) = 1; % 将细化线映射到边缘图像
figure
imshow(showimg),title('细化去支前');

% 细化线的交点、端点检测
% 端点坐标
endpox = zeros(1,lenbw4); % 不清楚大小，故设为lenbw4
endpoy = zeros(1,lenbw4);
% 交点坐标
intpox = zeros(1,lenbw4);
intpoy = zeros(1,lenbw4);
j = 1;k = 1;
for i = 1:lenbw4
    Rect = bw4((Ybw4(i)-1):(Ybw4(i)+1),(Xbw4(i)-1):(Xbw4(i)+1)); % 3X3矩形
    SumRect = sum(sum(Rect)) ; % 3X3对矩形求和
    if(SumRect == 2) % 端点条件
        endpox(j) = Xbw4(i);
        endpoy(j) = Ybw4(i);
        j = j + 1;
    elseif((SumRect >= 4)) % 交点条件
        intpox(k) = Xbw4(i);
        intpoy(k) = Ybw4(i);
        k = k + 1;
    end
end
% 原始端点和交点
indexend = find(endpox~=0); % 非零索引
indexint = find(intpox~=0);
% 交点、端点去重 后向比较
thr = 5; % 判断重复的阈值
for i = 1: length(indexend)-1 % 至少两个端点
    if(endpox(i)~=0) % 排除零点
        for j = i+1: length(indexend)
            d = sqrt((endpox(i)-endpox(j)).^2 + (endpoy(i)-endpoy(j)).^2);
            if(d<thr)
                endpox(j) = 0; % 距离太近的赋值0
                endpoy(j) = 0;
            end
        end
    end
end
% 可能没有交点或者一个交点
if(length(indexint)<=1)
else % 两个交点以上
    for i = 1:length(indexint)-1
        if(intpox(i)~=0) % 排除零点
            for j = i+1:length(indexint)
                d = sqrt((intpox(i)-intpox(j)).^2 + (intpoy(i)-intpoy(j)).^2);
                if(d<thr)
                    intpox(j) = 0;
                    intpoy(j) = 0;
                end
            end
        end
    end
end
indexend = find(endpox~=0); % 非零索引
indexint = find(intpox~=0);
lonend = length(indexend);
lonint = length(indexint);
% 将非零点集中到前端
for i = 1:lonint
    if(indexint(i)~=i)
        intpox(i) = intpox(indexint(i)) ;
        intpoy(i) = intpoy(indexint(i)) ;
        intpox(indexint(i)) = 0;
        intpoy(indexint(i)) = 0;
    end
end
intpox = intpox(1:length(indexint));
intpoy = intpoy(1:length(indexint));
for i = 1:lonend
    if(indexend(i)~=i)
        endpox(i) = intpox(indexend(i)) ;
        endpoy(i) = intpoy(indexend(i)) ;
        endpox(indexend(i)) = 0;
        endpoy(indexend(i)) = 0;
    end
end
endpox = endpox(1:length(indexend));
endpoy = endpoy(1:length(indexend));

% hold on
% plot(endpox,endpoy,'.','markersize',20);
% plot(intpox,intpoy,'*','markersize',10);


% 把检测到的交点背景化，使线段离散化
bw5 = bw4;
for i = 1:lonint
    bw5(intpoy(i)-1:intpoy(i)+1,intpox(i)-1:intpox(i)+1) = 0;
end
[Map,num] = bwlabel(bw5,8); % 标记
% figure
% showimg2 = bwedge;
% showimg2(bw5) = 1;
% imshow(bw4),title('分支离散化、求长度')

% 标记各交叉点、端点
% IntMark 数据格式
% 点编号
% X坐标
% Y坐标
% 交点周围三根线的编号
IntMark = zeros(6,lonint);
IntMark(2,:) = intpox;
IntMark(3,:) = intpoy;
for i = 1:lonint
    IntMark(1,i) = i;
    Rect = Map(intpoy(i)-2:intpoy(i)+2,intpox(i)-2:intpox(i)+2);
    temp= unique(Rect);
    IntMark(4:6,i) = temp(2:4,1);
end
% EndMark 数据格式
% 点编号
% X坐标
% Y坐标
% 端点所处的线段编号
EndMark = zeros(4,lonend);
EndMark(2,:) = endpox;
EndMark(3,:) = endpoy;
for i = 1:lonend
    EndMark(1,i) = i;
    Rect = Map(endpoy(i)-2:endpoy(i)+2,endpox(i)-2:endpox(i)+2);
    temp= unique(Rect);
    EndMark(4,i) = temp(2,1);
end
% 在图像中标记各线段标号
for i = 1:length(EndMark(1,:))
    text(EndMark(2,i),EndMark(3,i),num2str(EndMark(4,i)),'horiz','center','color','r','fontsize',25)
end

% 计算各分支长度并标记
% LineMark 数据格式 从大到小排序
% 线段编号
% 长度（点数）
LineMark = zeros(2,num);
maxcont = -inf ;
flag = 0; % 最长线段的编号
Cont = zeros(1,num);
for k = 1:num
    for i = 1:h
        for j = 1:w
            if(Map(i,j)==k)
                Cont(k) = Cont(k) + 1;
                LineMark(1,k) = k;
                LineMark(2,k) = Cont(k);
            end
        end
    end
    if(maxcont<Cont(k))
        maxcont = Cont(k);
        flag = k;
    end
end
for i = 1:length(LineMark(1,:))-1
    for j = i+1:length(LineMark(1,:))
        temp = zeros(2,1);
        if(LineMark(2,i) < LineMark(2,j))
            temp = LineMark(:,j);
            LineMark(:,j) = LineMark(:,i);
            LineMark(:,i) = temp;
        end
    end
end

% 只保留尾端最长支
endbw = bw5;
for i = 1:h
    for j = 1:w
        if(Map(i,j) ~= flag)
            endbw(i,j) = 0;
        end
    end
end
% 计算各分支对边缘图像的最小内切圆
% R 数据格式
% 线段编号
% 最小内切圆
R = zeros(2,num);
for n = 1:num
    minx = inf;
    for i = 1:h
        for j = 1:w
            if(Map(i,j) == n)             
                for k = 1:lenedge
                    d = sqrt((j-Xedge(k)).^2 + (i-Yedge(k)).^2);
                    if(minx > d)
                        minx = d;
                    end
                end
            end
        end
    end
    R(1,n) = n;
    R(2,n) = minx;
end
quzhibw = endbw; % copy尾端骨线
for n = 1:num
    % 1、不是尾端线段  2、最小内切圆半径>20（即不接触边缘）  3、线段长度大于25
    % 满足三个条件的线段保留
    if(R(1,n)~=flag && R(2,n)>20 && LineMark(2,find(LineMark(1,:)==n))>25)
        for i = 1:h
            for j = 1:w
                if(Map(i,j) == n)
                    quzhibw(i,j) = 1;
                end
            end
        end     
    end
end
indexquzhi = find(quzhibw ~= 0); % 非零点坐标索引
[rowquzhi,colquzhi] = ind2sub(size(quzhibw),indexquzhi);
Yquzhi = rowquzhi';
Xquzhi = colquzhi';
lenquzhi = length(Xquzhi);


% 对尾部图像处理
indexendbw = find(endbw ~= 0); % 非零点坐标索引
[rowend,colend] = ind2sub(size(endbw),indexendbw);
Yend = rowend';
Xend = colend';
lenend = length(Xend);
% 记录尾部两端点的属性
rec = zeros(1,2);
j = 1;
for i = 1:lonint
    if(any(IntMark(4:6,i)==flag))
        rec(j) = i;
        j = j + 1;
    end
end
weiendpo = zeros(2,2); % 两端点坐标
weiend = zeros(2,1); % 尾部起点（距离质心较远的端点）
j = 1;
de = zeros(1,2);
if(length(find(rec==0))==2) % 两端点
    Factor = 3.0; % 直径和长度的比例因子
    for i = 1:length(EndMark(1,:))
        if(EndMark(4,i)==flag)
            weiendpo(1,j) = EndMark(2,i) ;
            weiendpo(2,j) = EndMark(3,i) ;
            de(j) = sqrt((Xcenter-weiendpo(1,j)).^2 + (Ycenter-weiendpo(2,j)).^2);
            j = j + 1;
        end
    end
    Xran = abs(weiendpo(1,1)-weiendpo(1,2));
    Yran = abs(weiendpo(2,1)-weiendpo(2,2));
    if(de(1) > de(2))
        weiend(1) = weiendpo(1,1);
        weiend(2) = weiendpo(2,1);
    else
        weiend(1) = weiendpo(1,2);
        weiend(2) = weiendpo(2,2);
    end
elseif(length(find(rec==0))==1) % 一个端点，一个交叉点
    Factor = 3.2;
    weiendpo(1,1) = IntMark(2,rec(1)) ;
    weiendpo(2,1) = IntMark(3,rec(1)) ;
    for i = 1:length(EndMark(1,:))
        if(EndMark(4,i)==flag)
            weiendpo(1,2) = EndMark(2,i) ;
            weiendpo(2,2) = EndMark(3,i) ;
            weiend(1) = weiendpo(1,2);
            weiend(2) = weiendpo(2,2);
        end
    end
    Xran = abs(weiendpo(1,1)-weiendpo(1,2));
    Yran = abs(weiendpo(2,1)-weiendpo(2,2));
else % 两个交叉点，对尾部交叉点的两分支求平均
    Factor = 3.0;
    text(Xcenter,Ycenter,num2str(flag),'horiz','center','color','r','fontsize',25)
    weiendpo(1,1) = IntMark(2,rec(1)) ;
    weiendpo(2,1) = IntMark(3,rec(1)) ;
    weiendpo(1,2) = IntMark(2,rec(2)) ;
    weiendpo(2,2) = IntMark(3,rec(2)) ;
    Xran = abs(weiendpo(1,1)-weiendpo(1,2));
    Yran = abs(weiendpo(2,1)-weiendpo(2,2));
    d1 = sqrt((Xcenter-IntMark(2,rec(1))).^2 + (Ycenter-IntMark(3,rec(1))).^2);
    d2 = sqrt((Xcenter-IntMark(2,rec(2))).^2 + (Ycenter-IntMark(3,rec(2))).^2);
    sumx = 0;
    sumy = 0;
    if(d1 > d2)    
        for i = 1:length(EndMark(1,:))
             if(IntMark(4,rec(1))==EndMark(4,i) || IntMark(5,rec(1))==EndMark(4,i) || IntMark(6,rec(1))==EndMark(4,i))
                sumx = sumx + EndMark(2,i);
                sumy = sumy + EndMark(3,i);
            end
        end
    else
        for i = 1:length(EndMark(1,:))
            if(IntMark(4,rec(2))==EndMark(4,i) || IntMark(5,rec(2))==EndMark(4,i) || IntMark(6,rec(2))==EndMark(4,i))
                sumx = sumx + EndMark(2,i);
                sumy = sumy + EndMark(3,i);
            end
        end
    end
    weiend(1) = sumx/2;
    weiend(2) = sumy/2;
%     quzhibw = DrawLineImage(quzhibw,weiend(1),weiend(2),weiendpo(1,1),weiendpo(2,1));
end
figure
showquzhi = bwedge;
showquzhi(quzhibw) = 1;
imshow(showquzhi),title('细化去支后');
hold on
plot(weiendpo(1,1),weiendpo(2,1),'*','markersize',10);
plot(weiendpo(1,2),weiendpo(2,2),'*','markersize',10);
plot(weiend(1),weiend(2),'*','markersize',20);

xcon = (weiendpo(1,1)+weiendpo(1,2))/2;
ycon = (weiendpo(2,1)+weiendpo(2,2))/2;
if(Xran > Yran)
    yen = find(endbw(:,floor(xcon)));
    if(yen < ycon)
        Dre = 'up';
    else
        Dre = 'down';
    end
else
    xen = find(endbw(floor(ycon),:));
    if(xen < xcon)
        Dre = 'left';
    else
        Dre = 'right';
    end
end

% 求尾部骨线上的内切圆
Rr = zeros(1,lenend);
for o = 1:lenend
    minx = inf;
    for i = 1:lenedge      
        d = sqrt((Xend(o)-Xedge(i)).^2 + (Yend(o)-Yedge(i)).^2);
        if(minx > d)
            minx = d;
        end
    end
    Rr(o) = minx;
end

% 划分区间，统计内切圆半径数目
step = 1; % 区间长度
N = ceil(max(Rr)/step);
Reg = zeros(2,N);
for n = 1:N
    Reg(1,n) = step * n;
    for o = 1:lenend
        if(Rr(o) > step*n && Rr(o) < step*(n+1))
           Reg(2,n) = Reg(2,n)+1;
        end
    end
end

% % 显示
% figure('NumberTitle', 'off', 'Name', '尾部内切圆半径统计');
% plot(Reg(1,:),Reg(2,:),'*','markersize',10);
% xlabel('半径'),ylabel('次数'),title('尾部内切圆半径统计');

% 防止因为尾部骨线过长，使得统计出的半径偏向于头部
maxposition = zeros(1,2);
maxposition(1,1) = find(Reg(2,:) == max(Reg(2,:)));
Count = Reg(2,maxposition(1,1));
Reg(2,maxposition(1,1)) = 0;
maxposition(1,2) = find(Reg(2,:) == max(Reg(2,:)));
if(maxposition(1,1) - maxposition(1,2)>=2)
    Dmax = 2 * Reg(1,maxposition(1,2));
    Count = Reg(2,maxposition(1,2));
else
    Dmax = 2 * Reg(1,maxposition(1,1));
end

L1 = Factor * Dmax; % 尾翅长度


for alpha=0:pi/90:2*pi;
    x = weiend(1) + L1*cos(alpha);
    y = weiend(2) + L1*sin(alpha);
%     plot(x,y,'.')
    fx = floor(x);
    fy = floor(y);
    if(fx>1 && fx<w && fy>1 && fy<h)
        Rect = endbw(fy-1:fy+1,fx-1:fx+1);
        if(sum(sum(Rect))>=1)
            xc = fx;
            yc = fy;
        end
    end
end

L2 = 1/1.5 * L1; % 最后一关节长度

for alpha=0:pi/90:2*pi;
    x = xc + L2*cos(alpha);
    y = yc + L2*sin(alpha);
%     plot(x,y,'.')
    fx = floor(x);
    fy = floor(y);
    if(fx>1 && fx<w && fy>1 && fy<h)
        Rect = endbw(fy-1:fy+1,fx-1:fx+1);
        if(sum(sum(Rect))>=1)
            d = sqrt((fx-weiend(1)).^2 + (fy-weiend(2)).^2);
            if(d>L1)
                xc2 = fx;
                yc2 = fy;
            end
        end
    end
end

% 求整体骨线上最大半径和位置
Rr = zeros(1,lenquzhi);
for o = 1:lenquzhi
    minn = inf;
    for i = 1:lenedge      
        d = sqrt((Xquzhi(o)-Xedge(i)).^2 + (Yquzhi(o)-Yedge(i)).^2);
        if(minn > d)
            minn = d;
        end
    end
    Rr(o) = minn;
end
indexRr = find(Rr==max(Rr));
BanX = Xquzhi(indexRr);
BanY = Yquzhi(indexRr);
MaxR = Rr(indexRr);

figure
r = img(:,:,1);
g = img(:,:,2);
b = img(:,:,3);
r(quzhibw) = 255;
g(quzhibw) = 0;
b(quzhibw) = 0;
img(:,:,1) = r ;
img(:,:,2) = g ;
img(:,:,3) = b ;
imshow(img);
hold on
plot(xc,yc,'*','markersize',20,'color','r'); % 绘制尾翅位置
plot(xc2,yc2,'*','markersize',20,'color','g'); % 绘制最后一关节位置
plot(BanX,BanY,'*','markersize',20,'color','b'); % 绘制最大圆圆心
alpha=0:pi/20:2*pi;
x = BanX + MaxR*cos(alpha);
y = BanY + MaxR*sin(alpha);
plot(x,y,'.','color','b')
% %多项式拟合
% n=3;
% if(strcmp(Dre,'up')||strcmp(Dre,'down'))
%     A=polyfit(Xquzhi,Yquzhi,n);  %n是给定的多项式的次数，拟合出来的结果A是系数向量
%     y1=polyval(A,1:w);  %计算出拟合的y值
%     plot(1:w,y1,'r-');
% end
% if(strcmp(Dre,'left')||strcmp(Dre,'right'))
%     A=polyfit(Yquzhi,Xquzhi,n);  %n是给定的多项式的次数，拟合出来的结果A是系数向量
%     y1=polyval(A,1:h);  %计算出拟合的y值
%     plot(1:h,y1,'r-');
% end
% % 设出圆锥曲线方程
% F=@(p,x)p(1)*x(:,1).^2+p(2)*x(:,1).*x(:,2)+p(3)*x(:,2).^2+p(4);
% % 离散数据点
% x = zeros(length(Xquzhi),2);
% x(:,1) = Xquzhi';
% x(:,2) = Yquzhi';
% p0=[1 1 1 1];
% warning off
% % 拟合系数，最小二乘方法
% p=nlinfit(x,zeros(size(x,1),1),F,p0);
% plot(x(:,1),x(:,2),'ro');
% hold on;
% xmin=min(x(:,1));
% xmax=max(x(:,1));
% ymin=min(x(:,2));
% ymax=max(x(:,2));
% % 作图
% ezplot(@(x,y)F(p,[x,y]),[-1+xmin,1+xmax,-1+ymin,1+ymax]);
% title('曲线拟合');
% legend('样本点','拟合曲线')




