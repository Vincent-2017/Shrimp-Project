clear ;clc ;close all ;

img=imread('1.jpg');
[h,w,n] = size(img);
hsv=rgb2hsv(img); % 转hsv
s=hsv(:,:,2); % s分量
bw=1-im2bw(s,graythresh(s)); % 二值化
se=strel('disk',5);
bw2=imclose(bw,se); % 闭操作:通常消弥狭窄的间断和长细的鸿沟，消除小的空洞，并填补轮廓线中的断裂。
bw3=bwareaopen(bw2,500); % 删除二值图像BW中面积小于200的对象
figure(1)
subplot(221),imshow(bw),title('二值化');
subplot(222),imshow(bw2),title('闭操作');
subplot(223),imshow(bw3),title('去孤岛');
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
subplot(224),imshow(bw5),title('骨架去支');

hold on
alpha=0:pi/20:2*pi; 
R=100;         
x=R*cos(alpha); 
y=R*sin(alpha); 
plot(x,y,'-') 
fill(x,y,'r');   
%perim=bwperim(bw3); % 查找二值图像的边缘
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

