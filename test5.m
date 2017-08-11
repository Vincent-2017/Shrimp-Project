clear ;clc ;close all ;
% ��ֵ��
img=imread('1.jpg');
[h,w,n] = size(img);
hsv=rgb2hsv(img);  % תhsv
s=hsv(:,:,2);      % s����
bw=1-im2bw(s,graythresh(s)); % ��ֵ��
se=strel('disk',5);
bw2=imclose(bw,se);
bw3=bwareaopen(bw2,500);

% ���Ե
bbw3=double(bw3);
bwedge=edge(bbw3);
indexedge = find(bwedge ~= 0);
[rowedge,coledge] = ind2sub(size(bwedge),indexedge);
Yedge = rowedge';
Xedge = coledge';
lenedge = length(Yedge);

% ���Եͼ�������
Xsum = 0;
Ysum = 0;
for i = 1:lenedge
    Xsum = Xsum + Xedge(i);
    Ysum = Ysum + Yedge(i);
end
Xcenter = Xsum / lenedge ;
Ycenter = Ysum / lenedge ;

% ϸ��
bw4 = bwmorph(bw3,'thin',Inf);
indexbw4 = find(bw4 ~= 0); % �������������
[rowbw4,colbw4] = ind2sub(size(bw4),indexbw4);
Ybw4 = rowbw4';
Xbw4 = colbw4';
lenbw4 = length(Xbw4);

showimg = bwedge;
showimg(bw4) = 1; % ��ϸ����ӳ�䵽��Եͼ��
figure
imshow(showimg),title('���㡢�˵���');

% ϸ���ߵĽ��㡢�˵���
% �˵�����
endpox = zeros(1,lenbw4); % �������С������Ϊlenbw4
endpoy = zeros(1,lenbw4);
% ��������
intpox = zeros(1,lenbw4);
intpoy = zeros(1,lenbw4);
j = 1;k = 1;
for i = 1:lenbw4
    Rect = bw4((Ybw4(i)-1):(Ybw4(i)+1),(Xbw4(i)-1):(Xbw4(i)+1)); % 3X3����
    SumRect = sum(sum(Rect)) ; % 3X3�Ծ������
    if(SumRect == 2) % �˵�����
        endpox(j) = Xbw4(i);
        endpoy(j) = Ybw4(i);
        j = j + 1;
    elseif((SumRect >= 4)) % ����������������Ϊ (SumRect == 4)?
        intpox(k) = Xbw4(i);
        intpoy(k) = Ybw4(i);
        k = k + 1;
    end
end
indexend = find(endpox~=0); % ��������
indexint = find(intpox~=0);

% ���㡢�˵�ȥ��
thr = 5; %�жϾ������ֵ
for i = 1: length(indexend)-1 % ���������˵�
    if(endpox(i)~=0)
        for j = i+1: length(indexend)
            d = sqrt((endpox(i)-endpox(j)).^2 + (endpoy(i)-endpoy(j)).^2);
            if(d<thr)
                endpox(j) = 0; % ����̫���ĸ�ֵ0
                endpoy(j) = 0;
            end
        end
    end
end
% ����û�н������һ������
if(length(indexint)<=1)
else % ������������
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
indexend = find(endpox~=0); % ��������
indexint = find(intpox~=0);
% max = -inf ;
% for i = length(indexend)
%     d = sqrt((Xcenter-endpox(i)).^2 + (Ycenter-endpoy(i)).^2);
%     if(max < d)
%         max = d ;
%         maxx = endpox(i);
%         maxy = endpoy(i);
%     end
% end

hold on
plot(endpox(indexend),endpoy(indexend),'.','markersize',20);
plot(intpox(indexint),intpoy(indexint),'*','markersize',10);


% �Ѽ�⵽�Ľ��㱳����
for i = 1:length(indexint)
    bw4(intpoy(indexint(i))-1:intpoy(indexint(i))+1,intpox(indexint(i))-1:intpox(indexint(i))+1) = 0;
end
[Map,num] = bwlabel(bw4,8);
figure
showimg2 = bwedge;
showimg2(bw4) = 1;
imshow(showimg2),title('��֧��ɢ�����󳤶�')

% �������֧����
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

% ֻ����β���֧
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
imshow(showimg3),title('�������β��')

% ��β��ͼ����˵�
indexend = find(endbw ~= 0); % �������������
[rowend,colend] = ind2sub(size(endbw),indexend);
Yend = rowend';
Xend = colend';
lenend = length(Xend);
% �˵�����
weiendpox = zeros(1,2);
weiendpoy = zeros(1,2);
j = 1;
for i = 1:lenend
    Rect = endbw((Yend(i)-1):(Yend(i)+1),(Xend(i)-1):(Xend(i)+1)); % 3X3����
    SumRect = sum(sum(Rect)) ; % 3X3�Ծ������
    if(SumRect == 2) % �˵�����
        weiendpox(j) = Xend(i);
        weiendpoy(j) = Yend(i);
        j = j + 1;
    end
end
% hold on
% plot(weiendpox,weiendpoy,'*','markersize',10);
Xran = abs(weiendpox(1)-weiendpox(2));
Yran = abs(weiendpoy(1)-weiendpoy(2));
if(Xran>Yran)
    Dre = 'horizontal';
else
    Dre = 'vertical';
end

% ��β�����
max = -inf;
for i = 1:lenend
    d = sqrt((Xcenter-Xend(i)).^2 + (Ycenter-Yend(i)).^2);
    if(max < d)
        max = d ;
        xend = Xend(i);
        yend = Yend(i);
    end
end
% hold on
% plot(Xcenter,Ycenter,'.','markersize',40);
% plot(xend,yend,'.','markersize',20);
% line([Xcenter,xend],[Ycenter,yend])



Rr = zeros(1,lenend);
% ��β������Բ
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
% % ����β������Բ
% figure
% imshow(bwedge),title('����Բ');
% hold on
% alpha=0:pi/20:2*pi;
% for o = 1:lenend
%     x=Xend(o) + Rr(o)*cos(alpha);
%     y=Yend(o) + Rr(o)*sin(alpha);
%     plot(x,y,'.')
%     fill(x,y,'r');
% end

% ���ư뾶����
figure('NumberTitle', 'off', 'Name', '����������Բ�뾶');
if strcmp(Dre,'horizontal')
    o = 1:lenend;
    plot(Xend(o),Rr(o),'.')
    xlabel('X'),ylabel('R'),title('β������������Բ�뾶');
end
if strcmp(Dre,'vertical')
    o = 1:lenend;
    plot(Yend(o),Rr(o),'.')
    xlabel('Y'),ylabel('R'),title('β������������Բ�뾶');
end
