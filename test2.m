clear ;clc ;close all ;
% ��ֵ��
img=imread('3.jpg');
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
indexbw4 = find(bw4 ~= 0);
[rowbw4,colbw4] = ind2sub(size(bw4),indexbw4);
Ybw4 = rowbw4';
Xbw4 = colbw4';
lenbw4 = length(Xbw4);
max1 = -inf;
for i = 1:lenbw4
    d = sqrt((Xcenter-Xbw4(i)).^2 + (Ycenter-Ybw4(i)).^2);
    if(max1 < d)
        max1 = d ;
        lonx = Xbw4(i);
        lony = Ybw4(i);
    end
end
showimg = bwedge;
showimg(bw4) = 1;
figure
imshow(showimg),title('���㡢�˵���');

% ���㡢�˵���
endpox = zeros(1,lenbw4);
endpoy = zeros(1,lenbw4);
intpox = zeros(1,lenbw4);
intpoy = zeros(1,lenbw4);
endd = zeros(1,lenbw4);
intd = zeros(1,lenbw4);
j = 1;k = 1;
for i = 1:lenbw4
    Rect = bw4((Ybw4(i)-1):(Ybw4(i)+1),(Xbw4(i)-1):(Xbw4(i)+1));
    SumRect = sum(sum(Rect)) ;
    if(SumRect == 2)
        endpox(j) = Xbw4(i);
        endpoy(j) = Ybw4(i);
        endd(j) = 100;
        if(j>1)
            endd(j) = sqrt((endpox(j)-endpox(j-1)).^2 + (endpoy(j)-endpoy(j-1)).^2);
        end
        j = j + 1;
    elseif((SumRect == 4))
        intpox(k) = Xbw4(i);
        intpoy(k) = Ybw4(i);   
        intd(1) = 100;
        if(k>1)
            intd(k) = sqrt((intpox(k)-intpox(k-1)).^2 + (intpoy(k)-intpoy(k-1)).^2);
        end
        k = k + 1;
        
    end   
end
indexintd = find(intd>5);
indexendd = find(endd>5);
hold on
plot(endpox(indexendd),endpoy(indexendd),'.','markersize',15);
plot(intpox(indexintd),intpoy(indexintd),'*','markersize',10);
