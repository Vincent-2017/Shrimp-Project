% 原理：细化后的线由单点组成，根据点周围的状态不同来确定属性 
% 0 0 0
% 0 1 0 端点 和为2
% 0 0 1

% 0 0 0
% 1 1 0 一般点 和为3
% 0 0 1

% 0 1 0
% 1 1 0 交点 和为4（或者及以上？）
% 0 0 1

% 问题：交叉点附近可能会重复检测，通过判断检测到的点间的距离来排除同位置的交叉点（也对端点进行了同样处理，感觉不用也可以）

clear ;clc ;close all ;
% 二值化
img=imread('4.jpg');
[h,w,n] = size(img);
hsv=rgb2hsv(img);  % 转hsv
s=hsv(:,:,2);      % s分量
bw=1-im2bw(s,graythresh(s)); % 二值化
se=strel('disk',10);
bw2=imclose(bw,se);       % 闭操作:通常消弥狭窄的间断和长细的鸿沟，消除小的空洞，并填补轮廓线中的断裂。
bw3=bwareaopen(bw2,500);  % 删除二值图像BW中面积小于500的对象

% 求边缘
bbw3=double(bw3);
bwedge=edge(bbw3);

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

% 交点、端点检测
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
hold on
plot(endpox(indexend),endpoy(indexend),'.','markersize',20);
plot(intpox(indexint),intpoy(indexint),'*','markersize',10);