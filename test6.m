clear ;clc ;close all ;
% 二值化
img = imread('2.jpg');
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
% 可能没有交点或者一个交点
if(length(indexint)<=1)
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

hold on
plot(endpox,endpoy,'.','markersize',20);
plot(intpox,intpoy,'*','markersize',10);


% 把检测到的交点背景化
for i = 1:lonint
    bw4(intpoy(i)-1:intpoy(i)+1,intpox(i)-1:intpox(i)+1) = 0;
end
[Map,num] = bwlabel(bw4,8);
% figure
% showimg2 = bwedge;
% showimg2(bw4) = 1;
% imshow(showimg2),title('分支离散化、求长度')

% 标记各交叉点、端点
% 点编号
% X坐标
% Y坐标
% 附加属性
XIntMark = zeros(6,lonint);
XIntMark(2,:) = intpox;
XIntMark(3,:) = intpoy;
for i = 1:lonint
    XIntMark(1,i) = i;
    Rect = Map(intpoy(i)-2:intpoy(i)+2,intpox(i)-2:intpox(i)+2);
    temp= unique(Rect);
    XIntMark(4:6,i) = temp(2:4,1);
end
XEndMark = zeros(4,lonend);
XEndMark(2,:) = endpox;
XEndMark(3,:) = endpoy;
for i = 1:lonend
    XEndMark(1,i) = i;
    Rect = Map(endpoy(i)-2:endpoy(i)+2,endpox(i)-2:endpox(i)+2);
    temp= unique(Rect);
    XEndMark(4,i) = temp(2,1);
end

% 计算各分支长度
LineMark = zeros(2,num);
maxcont = -inf ;
flag = 0;
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
figure
showimg3 = bwedge;
showimg3(endbw) = 1;
imshow(showimg3),title('保留最长的尾部')

% 对尾部图像检测端点
indexendbw = find(endbw ~= 0); % 非零点坐标索引
[rowend,colend] = ind2sub(size(endbw),indexendbw);
Yend = rowend';
Xend = colend';
lenend = length(Xend);
% 端点坐标
weiendpo = zeros(2,2);
j = 1;
for i = 1:lenend
    Rect = endbw((Yend(i)-1):(Yend(i)+1),(Xend(i)-1):(Xend(i)+1)); % 3X3矩形
    SumRect = sum(sum(Rect)) ; % 3X3对矩形求和
    if(SumRect == 2) % 端点条件
        weiendpo(1,j) = Xend(i);
        weiendpo(2,j) = Yend(i);
        j = j + 1;
    end
end
Xran = abs(weiendpo(1,1)-weiendpo(1,2));
Yran = abs(weiendpo(2,1)-weiendpo(2,2));
if(Xran>Yran)
    Dre = 'horizontal';
else
    Dre = 'vertical';
end

rec = zeros(2,2);
j = 1;
for i = 1:lonint
    if(any(XIntMark(4:6,i)==flag))
        rec(1,j) = XIntMark(2,i);
        rec(2,j) = XIntMark(3,i);
        j = j + 1;
    end
end
Rr = zeros(1,lenend);

% 求尾部内切圆
for i = 1:lenend
    min = inf;
    for j = 1:lenedge
        d = sqrt((Xend(i)-Xedge(j)).^2 + (Yend(i)-Yedge(j)).^2);
        if(min > d)
            min = d ;
        end
    end
    Rr(i) = min;
end
% 绘制尾部内切圆
figure
imshow(bwedge),title('内切圆');
hold on
alpha=0:pi/20:2*pi;
for o = 1:lenend
    x=Xend(o) + Rr(o)*cos(alpha);
    y=Yend(o) + Rr(o)*sin(alpha);
    plot(x,y,'.')
    fill(x,y,'r');
end

% 绘制半径函数
figure('NumberTitle', 'off', 'Name', '骨线上内切圆半径');
if strcmp(Dre,'horizontal')
    o = 1:lenend;
    plot(Xend(o),Rr(o),'.')
    xlabel('X'),ylabel('R'),title('尾部骨线上内切圆半径');
end
if strcmp(Dre,'vertical')
    o = 1:lenend;
    plot(Yend(o),Rr(o),'.')
    xlabel('Y'),ylabel('R'),title('尾部骨线上内切圆半径');
end




