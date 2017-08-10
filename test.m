clear ;clc ;close all ;
% ��ֵ��
img=imread('12.jpg');
[h,w,n] = size(img);
hsv=rgb2hsv(img);  % תhsv
s=hsv(:,:,2);      % s����
bw=1-im2bw(s,graythresh(s)); % ��ֵ��
se=strel('disk',10);
bw2=imclose(bw,se);       % �ղ���:ͨ��������խ�ļ�Ϻͳ�ϸ�ĺ蹵������С�Ŀն�������������еĶ��ѡ�
bw3=bwareaopen(bw2,500);  % ɾ����ֵͼ��BW�����С��500�Ķ���

% ϸ��
bw4 = bwmorph(bw3,'thin',Inf);

% �õ�ϸ���ߵĶ˵�
bwpo=endpoints(bw4);
indexpo = find(bwpo ~= 0);
[rowpo,colpo] = ind2sub(size(bwpo),indexpo);
Ypo = rowpo';
Xpo = colpo';
lenpo = length(Xpo);

% ��˵��������
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

% ȷ��Ϻ������
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

% ϸ��ȥ֧
bw5 = xihuaquzhi(bw4,Dre);

% ��ʾ���
showbw4 = bwedge;
showbw5 = bwedge;
figure(1)
showbw4(bw4) = 1;
showbw5(bw5) = 1;
subplot(121),imshow(showbw4),title('ϸ��ȥ֧ǰ');
subplot(122),imshow(showbw5),title('ϸ��ȥ֧��');
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
% subplot(221),imshow(bw),title('��ֵ��');
% subplot(222),imshow(bw2),title('�ղ���');
% subplot(223),imshow(bw3),title('ȥ�µ�');
% subplot(224),imshow(bw4),title('�Ǽ�ϸ��');
% 
figure(3)
% bw3(bw5)=0;
imshow(bw3),title('����������Բ�仯');
% % c��Ӧˮƽλ�á�r��Ӧ��ֱλ��
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

figure('NumberTitle', 'off', 'Name', '����������Բ�뾶');
if(strcmp(Dre,'vd&up')||strcmp(Dre,'vd&down'))
    o = 1:lenbw;
    plot(o,Rr(o),'.')
    xlabel('x'),ylabel('y'),title('����������Բ�뾶');
end
if(strcmp(Dre,'hd&left')||strcmp(Dre,'hd&right'))
    o = 1:lenbw;
    plot(Rr(o),o,'.')
    xlabel('x'),ylabel('y'),title('����������Բ�뾶');
end

