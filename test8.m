clear ;clc ;close all ;
% ��ֵ��
img = imread('1.jpg');
[h,w,n] = size(img);
hsv = rgb2hsv(img);  % תhsv
s = hsv(:,:,2);      % s����
bw = ~im2bw(s,graythresh(s)); % ��ֵ��
se = strel('disk',10);
bw2 = imclose(bw,se);
bw3 = bwareaopen(bw2,500);

% ���Ե
bbw3 = double(bw3);
bwedge = edge(bbw3);
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
imshow(showimg),title('β��ʶ��');

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
lonend = length(indexend);
lonint = length(indexint);
% ������㼯�е�ǰ��
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


% �Ѽ�⵽�Ľ��㱳����
for i = 1:lonint
    bw4(intpoy(i)-1:intpoy(i)+1,intpox(i)-1:intpox(i)+1) = 0;
end
[Map,num] = bwlabel(bw4,8);
% figure
% showimg2 = bwedge;
% showimg2(bw4) = 1;
% imshow(showimg2),title('��֧��ɢ�����󳤶�')

% ��Ǹ�����㡢�˵�
% ����
% X����
% Y����
% ��������
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

% �������֧���Ȳ����
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
% figure
% showimg3 = bwedge;
% showimg3(endbw) = 1;
% imshow(showimg3),title('�������β��')

% ��β��ͼ����˵�
indexendbw = find(endbw ~= 0); % �������������
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
elseif(length(find(rec==0))==1) % һ���˵㣬һ����֧��
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
else % ���������
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

% ������Բ
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

% �������䣬ͳ�ư뾶��Ŀ
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

% % ��ʾ
% figure('NumberTitle', 'off', 'Name', 'β������Բ�뾶ͳ��');
% plot(Reg(1,:),Reg(2,:),'*','markersize',10);
% xlabel('�뾶'),ylabel('����'),title('β������Բ�뾶ͳ��');

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



