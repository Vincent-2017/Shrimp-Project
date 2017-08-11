% ԭ��ϸ��������ɵ�����ɣ����ݵ���Χ��״̬��ͬ��ȷ������ 
% 0 0 0
% 0 1 0 �˵� ��Ϊ2
% 0 0 1

% 0 0 0
% 1 1 0 һ��� ��Ϊ3
% 0 0 1

% 0 1 0
% 1 1 0 ���� ��Ϊ4�����߼����ϣ���
% 0 0 1

% ���⣺����㸽�����ܻ��ظ���⣬ͨ���жϼ�⵽�ĵ��ľ������ų�ͬλ�õĽ���㣨Ҳ�Զ˵������ͬ�������о�����Ҳ���ԣ�

clear ;clc ;close all ;
% ��ֵ��
img=imread('4.jpg');
[h,w,n] = size(img);
hsv=rgb2hsv(img);  % תhsv
s=hsv(:,:,2);      % s����
bw=1-im2bw(s,graythresh(s)); % ��ֵ��
se=strel('disk',10);
bw2=imclose(bw,se);       % �ղ���:ͨ��������խ�ļ�Ϻͳ�ϸ�ĺ蹵������С�Ŀն�������������еĶ��ѡ�
bw3=bwareaopen(bw2,500);  % ɾ����ֵͼ��BW�����С��500�Ķ���

% ���Ե
bbw3=double(bw3);
bwedge=edge(bbw3);

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

% ���㡢�˵���
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
hold on
plot(endpox(indexend),endpoy(indexend),'.','markersize',20);
plot(intpox(indexint),intpoy(indexint),'*','markersize',10);