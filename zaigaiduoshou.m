clear ;clc ;close all ;
% ��ֵ��
img = imread('8.jpg');
[h,w,n] = size(img);
hsv = rgb2hsv(img);  % תhsv
s = hsv(:,:,2);      % s����
bw = ~im2bw(s,graythresh(s)); % ��ֵ��
se = strel('disk',10);
bw2 = imclose(bw,se);
bw3 = bwareaopen(bw2,500); % ��ֵ�������ս��

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
Ycenter = Ysum / lenedge ; % ����

% matlabϸ��
bw4 = bwmorph(bw3,'thin',Inf); % ԭʼϸ��ͼ
indexbw4 = find(bw4 ~= 0); % �������������
[rowbw4,colbw4] = ind2sub(size(bw4),indexbw4);
Ybw4 = rowbw4';
Xbw4 = colbw4';
lenbw4 = length(Xbw4);

% ��ʾ
showimg = bwedge;
showimg(bw4) = 1; % ��ϸ����ӳ�䵽��Եͼ��
figure
imshow(showimg),title('ϸ��ȥ֧ǰ');

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
    elseif((SumRect >= 4)) % ��������
        intpox(k) = Xbw4(i);
        intpoy(k) = Ybw4(i);
        k = k + 1;
    end
end
% ԭʼ�˵�ͽ���
indexend = find(endpox~=0); % ��������
indexint = find(intpox~=0);
% ���㡢�˵�ȥ�� ����Ƚ�
thr = 5; % �ж��ظ�����ֵ
for i = 1: length(indexend)-1 % ���������˵�
    if(endpox(i)~=0) % �ų����
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
        if(intpox(i)~=0) % �ų����
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

% hold on
% plot(endpox,endpoy,'.','markersize',20);
% plot(intpox,intpoy,'*','markersize',10);


% �Ѽ�⵽�Ľ��㱳������ʹ�߶���ɢ��
bw5 = bw4;
for i = 1:lonint
    bw5(intpoy(i)-1:intpoy(i)+1,intpox(i)-1:intpox(i)+1) = 0;
end
[Map,num] = bwlabel(bw5,8); % ���
% figure
% showimg2 = bwedge;
% showimg2(bw5) = 1;
% imshow(bw4),title('��֧��ɢ�����󳤶�')

% ��Ǹ�����㡢�˵�
% IntMark ���ݸ�ʽ
% ����
% X����
% Y����
% ������Χ�����ߵı��
IntMark = zeros(6,lonint);
IntMark(2,:) = intpox;
IntMark(3,:) = intpoy;
for i = 1:lonint
    IntMark(1,i) = i;
    Rect = Map(intpoy(i)-2:intpoy(i)+2,intpox(i)-2:intpox(i)+2);
    temp= unique(Rect);
    IntMark(4:6,i) = temp(2:4,1);
end
% EndMark ���ݸ�ʽ
% ����
% X����
% Y����
% �˵��������߶α��
EndMark = zeros(4,lonend);
EndMark(2,:) = endpox;
EndMark(3,:) = endpoy;
for i = 1:lonend
    EndMark(1,i) = i;
    Rect = Map(endpoy(i)-2:endpoy(i)+2,endpox(i)-2:endpox(i)+2);
    temp= unique(Rect);
    EndMark(4,i) = temp(2,1);
end
% ��ͼ���б�Ǹ��߶α��
for i = 1:length(EndMark(1,:))
    text(EndMark(2,i),EndMark(3,i),num2str(EndMark(4,i)),'horiz','center','color','r','fontsize',25)
end

% �������֧���Ȳ����
% LineMark ���ݸ�ʽ �Ӵ�С����
% �߶α��
% ���ȣ�������
LineMark = zeros(2,num);
maxcont = -inf ;
flag = 0; % ��߶εı��
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
endbw = bw5;
for i = 1:h
    for j = 1:w
        if(Map(i,j) ~= flag)
            endbw(i,j) = 0;
        end
    end
end
% �������֧�Ա�Եͼ�����С����Բ
% R ���ݸ�ʽ
% �߶α��
% ��С����Բ
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
quzhibw = endbw; % copyβ�˹���
for n = 1:num
    % 1������β���߶�  2����С����Բ�뾶>20�������Ӵ���Ե��  3���߶γ��ȴ���25
    % ���������������߶α���
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
indexquzhi = find(quzhibw ~= 0); % �������������
[rowquzhi,colquzhi] = ind2sub(size(quzhibw),indexquzhi);
Yquzhi = rowquzhi';
Xquzhi = colquzhi';
lenquzhi = length(Xquzhi);


% ��β��ͼ����
indexendbw = find(endbw ~= 0); % �������������
[rowend,colend] = ind2sub(size(endbw),indexendbw);
Yend = rowend';
Xend = colend';
lenend = length(Xend);
% ��¼β�����˵������
rec = zeros(1,2);
j = 1;
for i = 1:lonint
    if(any(IntMark(4:6,i)==flag))
        rec(j) = i;
        j = j + 1;
    end
end
weiendpo = zeros(2,2); % ���˵�����
weiend = zeros(2,1); % β����㣨�������Ľ�Զ�Ķ˵㣩
j = 1;
de = zeros(1,2);
if(length(find(rec==0))==2) % ���˵�
    Factor = 3.0; % ֱ���ͳ��ȵı�������
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
elseif(length(find(rec==0))==1) % һ���˵㣬һ�������
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
else % ��������㣬��β������������֧��ƽ��
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
imshow(showquzhi),title('ϸ��ȥ֧��');
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

% ��β�������ϵ�����Բ
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

% �������䣬ͳ������Բ�뾶��Ŀ
step = 1; % ���䳤��
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

% ��ֹ��Ϊβ�����߹�����ʹ��ͳ�Ƴ��İ뾶ƫ����ͷ��
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

L1 = Factor * Dmax; % β�᳤��


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

L2 = 1/1.5 * L1; % ���һ�ؽڳ���

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

% ��������������뾶��λ��
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
plot(xc,yc,'*','markersize',20,'color','r'); % ����β��λ��
plot(xc2,yc2,'*','markersize',20,'color','g'); % �������һ�ؽ�λ��
plot(BanX,BanY,'*','markersize',20,'color','b'); % �������ԲԲ��
alpha=0:pi/20:2*pi;
x = BanX + MaxR*cos(alpha);
y = BanY + MaxR*sin(alpha);
plot(x,y,'.','color','b')
% %����ʽ���
% n=3;
% if(strcmp(Dre,'up')||strcmp(Dre,'down'))
%     A=polyfit(Xquzhi,Yquzhi,n);  %n�Ǹ����Ķ���ʽ�Ĵ�������ϳ����Ľ��A��ϵ������
%     y1=polyval(A,1:w);  %�������ϵ�yֵ
%     plot(1:w,y1,'r-');
% end
% if(strcmp(Dre,'left')||strcmp(Dre,'right'))
%     A=polyfit(Yquzhi,Xquzhi,n);  %n�Ǹ����Ķ���ʽ�Ĵ�������ϳ����Ľ��A��ϵ������
%     y1=polyval(A,1:h);  %�������ϵ�yֵ
%     plot(1:h,y1,'r-');
% end
% % ���Բ׶���߷���
% F=@(p,x)p(1)*x(:,1).^2+p(2)*x(:,1).*x(:,2)+p(3)*x(:,2).^2+p(4);
% % ��ɢ���ݵ�
% x = zeros(length(Xquzhi),2);
% x(:,1) = Xquzhi';
% x(:,2) = Yquzhi';
% p0=[1 1 1 1];
% warning off
% % ���ϵ������С���˷���
% p=nlinfit(x,zeros(size(x,1),1),F,p0);
% plot(x(:,1),x(:,2),'ro');
% hold on;
% xmin=min(x(:,1));
% xmax=max(x(:,1));
% ymin=min(x(:,2));
% ymax=max(x(:,2));
% % ��ͼ
% ezplot(@(x,y)F(p,[x,y]),[-1+xmin,1+xmax,-1+ymin,1+ymax]);
% title('�������');
% legend('������','�������')




