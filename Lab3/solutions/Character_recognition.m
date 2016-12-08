%%%% Character recognition 
clc, clear all; close all ; 
Path = 'images/';

img1 = mat2gray(im2double(rgb2gray(imread(fullfile(Path, 'a.png')))));
%figure; imhist(img1)
img2 = mat2gray(im2double(imread(fullfile(Path, 'text.png'))));
figure; imhist(img2)

mask1 = im2bw(img1, graythresh(img1));
mask2 = im2bw(img2, graythresh(img2));
mask1 = ~mask1; 

SE = strel('square', 2); 
mask2 = imdilate(mask2, SE); 
mask2 = imopen(mask2, SE); 

figure; subplot(121); imshow(mask1); subplot(122); imshow(mask2)
crr = xcorr2(im2double(mask2), im2double(mask1));


[ssr,snd] = sort(crr(:), 'descend');
ssr = ssr(1:20);
snd = snd(1:20);

[ij,ji] = ind2sub(size(crr),snd);

img3 = img2; 
figure ; imshow(img2); 
hold on ; 
for d = 1 : length(snd)
    plot(ji(d), ij(d), 'r*', 'MarkerSize', 10); 
    img3(ij(d):-1:ij(d)-size(mask1,1)+1,ji(d):-1:ji(d)-size(mask1,2)+1) = rot90(mask1,2);
end
hold off; 

figure; imshow(img3);

