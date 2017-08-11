clear ;clc ;close all ;

img=imread('1.jpg');
[h,w,n] = size(img);
hsv=rgb2hsv(img); % תhsv
s=hsv(:,:,2); % s����
bw=1-im2bw(s,graythresh(s)); % ��ֵ��
se=strel('disk',5);
bw2=imclose(bw,se); % �ղ���:ͨ��������խ�ļ�Ϻͳ�ϸ�ĺ蹵������С�Ŀն�������������еĶ��ѡ�
bw3=bwareaopen(bw2,500); % ɾ����ֵͼ��BW�����С��200�Ķ���
figure(1)
subplot(221),imshow(bw),title('��ֵ��');
subplot(222),imshow(bw2),title('�ղ���');
subplot(223),imshow(bw3),title('ȥ�µ�');
bw4 = bwmorph(bw3,'thin',Inf);
bw5 = bw4;
for j = 1:w
    if sum(bw4(:,j))>1
        for i = 1:h
            if bw4(i,j) == 1
                bw4(:,j) = 0;
                bw4(i,j) = 1;
                bw5(:,j) = bw4(:,j);
            end
        end     
    end
end
[r, c] = find(bw5 == 1);
subplot(224),imshow(bw5),title('�Ǽ�ȥ֧');

hold on
alpha=0:pi/20:2*pi; 
R=100;         
x=R*cos(alpha); 
y=R*sin(alpha); 
plot(x,y,'-') 
fill(x,y,'r');   
%perim=bwperim(bw3); % ���Ҷ�ֵͼ��ı�Ե
r=img(:,:,1);
g=img(:,:,2);
b=img(:,:,3);
r(bw5)=0;
g(bw5)=255;
b(bw5)=0;
img(:,:,1)=r;
img(:,:,2)=g;
img(:,:,3)=b;
figure
imshow(img)


        
% bw5 = imclearborder(bw2,8);
% imshow(bw5);

