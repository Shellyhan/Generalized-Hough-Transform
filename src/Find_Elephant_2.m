
clear;
close all;
%Read image:
img = imread('animals2.jpg');
temp = imread('template_elephant.png');
img_origin = img;

%get image gray scale:
img = rgb2gray(img);
temp = rgb2gray(temp);

%image processing for good edges:
img = medfilt2(img);
img = edge(img,'sobel',0.049);
temp = edge(temp,'sobel',0.045);

%reference zero point in template:
refPointx = round(size(temp,1)/2);
refPointy = round(size(temp,2)/2);

%-----------------------------
%Generate the R-table

%get template edge point:
[x,y] = find(temp > 0);
maxP = size(x,1);
%get the gradient of tamplate:
grad = Gradient(temp);
maxA = 180;

%Rtable with rotation and scale:
rtable = zeros(2*maxA, maxP, 2);
binCount = zeros(2*maxA,1);

for i=1:1:maxP
    k = grad(x(i), y(i)) + 180;
    binCount(k) = binCount(k) + 1;
    h = binCount(k);
%get dx and dy with scaling = 1.5:
    delta_x = 1.5*(x(i) - refPointx);
    delta_y = 1.5*(y(i) - refPointy);
%rotation:
    ang = -15*pi()/180;
    rtable(k, h, 1) = round(cos(ang)*delta_x - sin(ang)*delta_y);
    rtable(k, h, 2) = round(sin(ang)*delta_x + cos(ang)*delta_y);
end;

%-----------------------------
%Accumulator:

%get the image edge points
[a,b] = find(img > 0);
maxP_img =  size(a,1);

%gradient of image:
img_grad = Gradient(img);
size_img = size(img);

%generate accumulator:
count = zeros(size_img);
for i=1:1:maxP_img
    %the gradient angle:
    h = img_grad(a(i), b(i)) + 180;
    %count votes:
    for j = 1:1:binCount(h)
        c = a(i) - rtable(h, j, 1);
        d = b(i) - rtable(h, j, 2);
        if (c>0) && (c<size_img(1)) && (d>0) && (d<size_img(2))
            count(c, d) = count(c, d)+1;
        end;
    end;
end;
%-----------------------------

%Display matches:
subplot(221),imshow(img),title('Edges');

%find local maxima:
max_value = max(max(count));
[p,q] = find(count == max_value);

%the votes:
count1 = mat2gray(count);
subplot(222),imshow(count1),title('Accumulator Matrix');
hold on;
Circle(q, p, 0, 5);

subplot(223),imshow(img_origin),title('Result');
hold on;
plot(q,p,'g+', 'MarkerSize', 2);
hold on;
Circle(q, p, refPointy*1.5, 3);


function Circle(centery, centerx, reference, r)
radius = reference + r;
angle = 0:0.01:2*pi; 
d_x = radius*cos(angle);
d_y = radius*sin(angle);
plot(centery+d_y, centerx+d_x, 'g');
end




function [result] = Gradient(input)
    dy=imfilter(double(input),[1; -1],'same');
    dx=imfilter(double(input),[1  -1],'same');
    result = atan2(dy,dx)*180/pi();
end
