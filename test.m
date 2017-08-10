clear ;clc ;close all ;
% 二值化
img=imread('12.jpg');
[h,w,n] = size(img);
hsv=rgb2hsv(img);  % 转hsv
s=hsv(:,:,2);      % s分量
bw=1-im2bw(s,graythresh(s)); % 二值化
se=strel('disk',10);
bw2=imclose(bw,se);       % 闭操作:通常消弥狭窄的间断和长细的鸿沟，消除小的空洞，并填补轮廓线中的断裂。
bw3=bwareaopen(bw2,500);  % 删除二值图像BW中面积小于500的对象

% 细化
bw4 = bwmorph(bw3,'thin',Inf);

% 得到细化线的端点
bwpo=endpoints(bw4);
indexpo = find(bwpo ~= 0);
[rowpo,colpo] = ind2sub(size(bwpo),indexpo);
Ypo = rowpo';
Xpo = colpo';
lenpo = length(Xpo);

% 求端点间的最长连线
max = -inf;
for i = 1:lenpo
    for j = 1:lenpo
        d = sqrt((Xpo(i)-Xpo(j)).^2 + (Ypo(i)-Ypo(j)).^2);
        if(max < d)
            max = d ;
            x1 = Xpo(i);
            y1 = Ypo(i);
            x2 = Xpo(j);
            y2 = Ypo(j);
        end
    end
end
Xran = abs(x1 - x2);
Yran = abs(y1 - y2);
kmax = -1/((y2 - y1)/(x2 - x1));

% 求边缘
bbw3=double(bw3);
bwedge=edge(bbw3);
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

% 确定虾背方向
if(Xran > Yran)
    if(y1 > y2)
        Ycon = y1 - abs(y2 - y1)/abs(x2 - x1)*abs(x1 - Xcenter);
    else 
        Ycon = y1 + abs(y2 - y1)/abs(x2 - x1)*abs(x1 - Xcenter);
    end
    if(Ycon > Ycenter)
        Dre = 'vd&up';
    else
        Dre = 'vd&down';
    end
elseif(Xran < Yran)
    if(x1 > x2)
        Xcon = x1 - abs(x2 - x1)/abs(y2 - y1)*abs(y1 - Ycenter);
    else
        Xcon = x1 + abs(x2 - x1)/abs(y2 - y1)*abs(y1 - Ycenter);
    end
    if(Xcon > Xcenter)
        Dre = 'hd&left';
    else
        Dre = 'hd&right';
    end   
end

% 细化去支
bw5 = xihuaquzhi(bw4,Dre);

% 显示结果
showbw4 = bwedge;
showbw5 = bwedge;
figure(1)
showbw4(bw4) = 1;
showbw5(bw5) = 1;
subplot(121),imshow(showbw4),title('细化去支前');
subplot(122),imshow(showbw5),title('细化去支后');
hold on
plot(Xcenter,Ycenter,'.');
line([x1,x2],[y1,y2]);

% if(strcmp(Dre,'vd&up')||strcmp(Dre,'vd&down'))
%     line([Xcenter,Xcenter],[Ycenter,Ycon]);
% end
% if(strcmp(Dre,'hd&left')||strcmp(Dre,'hd&right'))
%      line([Xcenter,Xcon],[Ycenter,Ycenter]);
% end


% figure(2)
% subplot(221),imshow(bw),title('二值化');
% subplot(222),imshow(bw2),title('闭操作');
% subplot(223),imshow(bw3),title('去孤岛');
% subplot(224),imshow(bw4),title('骨架细化');
% 
figure(3)
% bw3(bw5)=0;
imshow(bw3),title('骨线上内切圆变化');
% % c对应水平位置、r对应垂直位置
indexbw = find(bw5 ~= 0);
[rowbw,colbw] = ind2sub(size(bw5),indexbw);
Ybw = rowbw';
Xbw = colbw';
lenbw = length(Xbw);
Rr = zeros(1,lenbw);


for o = 1:lenbw
    min = inf;
    for i = 1:lenedge
        d = sqrt((Xbw(o)-Xedge(i)).^2 + (Ybw(o)-Yedge(i)).^2);
        if(min > d)
            min = d;
        end
    end
    Rr(o) = min;
end
hold on
alpha=0:pi/20:2*pi; 
for o = 1:lenbw
    x=Xbw(o) + Rr(o)*cos(alpha);
    y=Ybw(o) + Rr(o)*sin(alpha);
    plot(x,y,'.')
    fill(x,y,'r');
end

figure('NumberTitle', 'off', 'Name', '骨线上内切圆半径');
if(strcmp(Dre,'vd&up')||strcmp(Dre,'vd&down'))
    o = 1:lenbw;
    plot(o,Rr(o),'.')
    xlabel('x'),ylabel('y'),title('骨线上内切圆半径');
end
if(strcmp(Dre,'hd&left')||strcmp(Dre,'hd&right'))
    o = 1:lenbw;
    plot(Rr(o),o,'.')
    xlabel('x'),ylabel('y'),title('骨线上内切圆半径');
end

