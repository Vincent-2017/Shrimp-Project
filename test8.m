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

showimg = bwedge;
showimg(bw4) = 1; % 将细化线映射到边缘图像
figure
imshow(showimg),title('尾翅识别');

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

% 标记各交叉点、端点
% 点编号
% X坐标
% Y坐标
% 附加属性
IntMark = zeros(6,lonint);
IntMark(2,:) = intpox;
IntMark(3,:) = intpoy;
for i = 1:lonint
    IntMark(1,i) = i;
    Rect = Map(intpoy(i)-2:intpoy(i)+2,intpox(i)-2:intpox(i)+2);
    temp= unique(Rect);
    IntMark(4:6,i) = temp(2:4,1);
end
EndMark = zeros(4,lonend);
EndMark(2,:) = endpox;
EndMark(3,:) = endpoy;
for i = 1:lonend
    EndMark(1,i) = i;
    Rect = Map(endpoy(i)-2:endpoy(i)+2,endpox(i)-2:endpox(i)+2);
    temp= unique(Rect);
    EndMark(4,i) = temp(2,1);
end
for i = 1:length(EndMark(1,:))
    text(EndMark(2,i),EndMark(3,i),num2str(EndMark(4,i)),'horiz','center','color','r','fontsize',25)
end

% 计算各分支长度并标记
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

% 对尾部图像检测端点
indexendbw = find(endbw ~= 0); % 非零点坐标索引
[rowend,colend] = ind2sub(size(endbw),indexendbw);
Yend = rowend';
Xend = colend';
lenend = length(Xend);

rec = zeros(1,2);
j = 1;
for i = 1:lonint
    if(any(IntMark(4:6,i)==flag))
        rec(j) = i;
        j = j + 1;
    end
end
weiendpo = zeros(2,2);
weiend = zeros(2,1);
j = 1;
de = zeros(1,2);
if(length(find(rec==0))==2)
    Factor = 3.0;
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
elseif(length(find(rec==0))==1) % 一个端点，一个分支点
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
else % 两个交叉点
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
    linepo = zeros(4,2);
    j = 1 ;
    if(d1 > d2)    
        for i = 1:length(EndMark(1,:))
             if(IntMark(4,rec(1))==EndMark(4,i) || IntMark(5,rec(1))==EndMark(4,i) || IntMark(6,rec(1))==EndMark(4,i))
                 linepo(1,j) = EndMark(4,i);
                 linepo(2,j) = EndMark(2,i);
                 linepo(3,j) = EndMark(3,i);
                 j = j + 1;
            end
        end
    else
        for i = 1:length(EndMark(1,:))
            if(IntMark(4,rec(2))==EndMark(4,i) || IntMark(5,rec(2))==EndMark(4,i) || IntMark(6,rec(2))==EndMark(4,i))
                 linepo(1,j) = EndMark(4,i);
                 linepo(2,j) = EndMark(2,i);
                 linepo(3,j) = EndMark(3,i);
                 j = j + 1;
            end
        end
    end
    for i = 1:num
        if(LineMark(1,i)==linepo(1,1))
            linepo(4,1) = LineMark(2,i);
        end
        if(LineMark(1,i)==linepo(1,2))
            linepo(4,2) = LineMark(2,i);
        end
    end
    if(linepo(4,2)>linepo(4,1))
        weiend(1) = linepo(2,2);
        weiend(2) = linepo(3,2);
    else
        weiend(1) = linepo(2,2);
        weiend(2) = linepo(3,2);
    end
        
end
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

% 求内切圆
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
step = 1;
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
Count = Reg(2,maxposition(1,1));
Reg(2,maxposition(1,1)) = 0;
maxposition(1,2) = find(Reg(2,:) == max(Reg(2,:)));
if(maxposition(1,1) - maxposition(1,2)>=2)
    Dmax = 2 * Reg(1,maxposition(1,2));
    Count = Reg(2,maxposition(1,2));
else
    Dmax = 2 * Reg(1,maxposition(1,1));
end

L = Factor * Dmax;

hold on
alpha=0:pi/90:2*pi;
x = weiend(1) + L*cos(alpha);
y = weiend(2) + L*sin(alpha);
plot(x,y,'.')



