
clear;
close all;
%Read image:
img = imread('animals2.jpg');
temp = imread('template_bear.png');
img_origin = img;

%modify image format:
img = rgb2gray(img);
temp = rgb2gray(temp);

%get image edge:
img = medfilt2(img);
img = medfilt2(img);
img = medfilt2(img);
img = medfilt2(img);
img = edge(img,'sobel',0.038);
temp = edge(temp,'canny');


%reference zero point:
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


%Rtable:
rtable = zeros(2*maxA, maxP, 2);
rtable_rotate = zeros(2*maxA, maxP, 2);
binCount = zeros(2*maxA,1);

%Rtable rotation:
for i=1:1:maxP
    
    k = grad(x(i), y(i)) + 180;
    binCount(k) = binCount(k) + 1;
    h = binCount(k);

    delta_x = 1.5*(x(i) - refPointx);
    delta_y = 1.5*(y(i) - refPointy);
%rotation-1:
    ang = -40*pi()/180;
    rtable(k, h, 1) = round(cos(ang)*delta_x - sin(ang)*delta_y);
    rtable(k, h, 2) = round(sin(ang)*delta_x + cos(ang)*delta_y);
%rotation-2:
    delta_x = 1.4*(x(i) - refPointx);
    delta_y = 1.4*(y(i) - refPointy);
    ang = -105*pi()/180;
    rtable_rotate(k, h, 1) = round(cos(ang)*delta_x - sin(ang)*delta_y);
    rtable_rotate(k, h, 2) = round(sin(ang)*delta_x + cos(ang)*delta_y);

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

%for rotation 1 match:
for i=1:1:maxP_img
    %the gradient angle:
    h = img_grad(a(i), b(i)) + 180;
    
    for j = 1:1:binCount(h)
        c = a(i) - rtable(h, j, 1);
        d = b(i) - rtable(h, j, 2);
        if (c>0) && (c<size_img(1)) && (d>0) && (d<size_img(2))
            count(c, d) = count(c, d)+1;
        end;
    end;
end;

%for rotation 2 match:
count_2 = zeros(size_img);

for i=1:1:maxP_img
    %the gradient angle:
    h = img_grad(a(i), b(i)) + 180;
    
    for j = 1:1:binCount(h)
        c = a(i) - rtable_rotate(h, j, 1);
        d = b(i) - rtable_rotate(h, j, 2);
        if (c>0) && (c<size_img(1)) && (d>0) && (d<size_img(2))
            count_2(c, d) = count_2(c, d)+1;
        end;
    end;
end;


%-----------------------------

%Display matches:
subplot(221),imshow(img),title('Edges');

%find local maxima:
max_value = max(max(count));
[p,q] = find(count == max_value);

max_2 = max(max(count_2));
[p1,q1] = find(count_2 == max_2);

%the votes:
count1 = mat2gray(count);
str = sprintf('Accumulator Matrix for object 1:(%d, %d)',q,p);
subplot(222),imshow(count1),title(str);
hold on;
Circle(q, p, 0, 5);

count2 = mat2gray(count_2);
str2 = sprintf('Accumulator Matrix for object 2:(%d, %d)',q1,p1);
subplot(223),imshow(count2),title(str2);
hold on;
Circle(q1, p1, 0, 5);

subplot(224),imshow(img_origin),title('Result');
hold on;
plot(q,p,'g+', 'MarkerSize', 2);
hold on;
plot(q1,p1,'g+', 'MarkerSize', 2);
hold on;
Circle(q, p, round(refPointy*1.5), 3);
hold on;
Circle(q1, p1, round(refPointy*1.5), 3);

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


