clear ;clc ;close all ;
% 二值化
img = imread('1.jpg');
[h,w,n] = size(img);
hsv = rgb2hsv(img);  % 转hsv
s = hsv(:,:,2);      % s分量
bw = ~im2bw(s,graythresh(s)); % 二值化
se = strel('disk',10);
bw2 = imclose(bw,se);
bw3 = bwareaopen(bw2,500);

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
Ycenter = Ysum / lenedge ;

% 细化
bw4 = bwmorph(bw3,'thin',Inf);
indexbw4 = find(bw4 ~= 0); % 非零点坐标索引
[rowbw4,colbw4] = ind2sub(size(bw4),indexbw4);
Ybw4 = rowbw4';
Xbw4 = colbw4';
lenbw4 = length(Xbw4);

% 显示
showimg = bwedge;
showimg(bw4) = 1; % 将细化线映射到边缘图像
figure
imshow(showimg),title('交点、端点检测');

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
    elseif((SumRect >= 4)) % 交点条件，或者设为 (SumRect == 4)?
        intpox(k) = Xbw4(i);
        intpoy(k) = Ybw4(i);
        k = k + 1;
    end
end
indexend = find(endpox~=0); % 非零索引
indexint = find(intpox~=0);
% 交点、端点去重
thr = 5; %判断距离的阈值
for i = 1: length(indexend)-1 % 至少两个端点
    if(endpox(i)~=0)
        for j = i+1: length(indexend)
            d = sqrt((endpox(i)-endpox(j)).^2 + (endpoy(i)-endpoy(j)).^2);
            if(d<thr)
                endpox(j) = 0; % 距离太近的赋值0
                endpoy(j) = 0;
            end
        end
    end
end
if(length(indexint)<=1) % 可能没有交点或者一个交点
else % 两个交点以上
    for i = 1:length(indexint)-1
        if(intpox(i)~=0)
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


% 把检测到的交点背景化
for i = 1:lonint
    bw4(intpoy(i)-1:intpoy(i)+1,intpox(i)-1:intpox(i)+1) = 0;
end
[Map,num] = bwlabel(bw4,8);
% figure
% showimg2 = bwedge;
% showimg2(bw4) = 1;
% imshow(showimg2),title('分支离散化、求长度')

% 计算各分支长度
% 最长分支标记为flag
maxcont = -inf ;
flag = 0;
Cont = zeros(1,num);
for k = 1:num
    for i = 1:h
        for j = 1:w
            if(Map(i,j)==k)
                Cont(k) = Cont(k) + 1;
            end
        end
    end
    if(maxcont<Cont(k))
        maxcont = Cont(k);
        flag = k;
    end
end

% 只保留尾端最长支
endbw = bw4;
cont = 0;
for i = 1:h
    for j = 1:w
        if(Map(i,j) ~= flag)
            endbw(i,j) = 0;
        end
    end
end
% figure
% showimg3 = bwedge;
% showimg3(endbw) = 1;
% imshow(showimg3),title('保留最长的尾部')

% 检测尾部端点
indexendbw = find(endbw ~= 0); % 非零点坐标索引
[rowend,colend] = ind2sub(size(endbw),indexendbw);
Yend = rowend';
Xend = colend';
lenend = length(Xend);
bwendpo = zeros(2,2);
j = 1;
for i = 1:lenend
    Rect = endbw((Yend(i)-1):(Yend(i)+1),(Xend(i)-1):(Xend(i)+1)); % 3X3矩形
    SumRect = sum(sum(Rect)) ; % 3X3对矩形求和
    if(SumRect == 2) % 端点条件
        bwendpo(1,j) = Xend(i);
        bwendpo(2,j) = Yend(i);
        j = j + 1;
    end
end
hold on
plot(bwendpo(1,:),bwendpo(2,:),'*','markersize',10);
% 确定尾部起点
weiend = zeros(2,1);
d1 = sqrt((Xcenter-bwendpo(1,1)).^2 + (Ycenter-bwendpo(2,1)).^2);
d2 = sqrt((Xcenter-bwendpo(1,2)).^2 + (Ycenter-bwendpo(2,2)).^2);
if(d1 > d2)
    weiend = bwendpo(:,1);
else
    weiend = bwendpo(:,2);
end

% 求尾部内切圆
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

% 划分区间，统计半径数目
step = 2;
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

maxposition = zeros(1,2);
maxposition(1,1) = find(Reg(2,:) == max(Reg(2,:)));
Reg(2,maxposition(1,1)) = 0;
maxposition(1,2) = find(Reg(2,:) == max(Reg(2,:)));
if(maxposition(1,1) - maxposition(1,2)>=2)
    Dmax = 2 * Reg(1,maxposition(1,2));
else
    Dmax = 2 * Reg(1,maxposition(1,1));
end

L = 3.2 * Dmax;

hold on
alpha=0:pi/90:2*pi;
x = weiend(1,1) + L*cos(alpha);
y = weiend(2,1) + L*sin(alpha);
plot(x,y,'.')
