rgbI = imread('Actin_5.jpg'); % Read the image file

grayImage = rgb2gray(rgbI); % Convert the colour image to a grayscale image
figure; imshow(grayImage); % Display grayscale image
figure; imhist(grayImage); % Display histogram of grayscale image

rImage = rgbI(:,:,1); % Separate the red component of the RGB image
figure; imshow(rImage); % Display image
figure; imhist(rImage); % Display histogram of image

gImage = rgbI(:,:,2); % Separate the green component of the RGB image
bImage = rgbI(:,:,3); % Separate the blue component of the RGB image

% Plot each image and it's histogram
figure; 
    subplot(2,4,1); imshow(grayImage); title('Grayscale');
    subplot(2,4,5); imhist(grayImage); [grayCounts,grayBinLocs] = imhist(grayImage); 
    subplot(2,4,2); imshow(rImage); title('Red');
    subplot(2,4,6); imhist(rImage); [rCounts,rBinLocs] = imhist(rImage); 
    subplot(2,4,3); imshow(gImage); title('Green');
    subplot(2,4,7); imhist(gImage); [gCounts,gBinLocs] = imhist(gImage); 
    subplot(2,4,4); imshow(bImage); title('Blue');
    subplot(2,4,8); imhist(bImage); [bCounts,bBinLocs] = imhist(bImage); 

J = imadjust(grayImage); % perform contrast stretching
figure; imshow(J); % Display image
    
otsuLevel = graythresh(bImage); % Get Otsu threshold for blue component
bwImage = im2bw(bImage,otsuLevel); % Convert blue image to black and white
figure; imshow(bwImage); % Display black and white image
    
se = strel('disk',5); % Create a disk shaped structuring element with 5 pixel radius
temp = imerode(bwImage,se); % Perform morphological erosion
nucleiImage = imclose(temp,se); % Perform morphological closing
figure; imshow(nucleiImage); % Display image

statsArea = regionprops(nucleiImage,'Area');
area = struct2array(statsArea);
figure; hist(area);

statsCentroid = regionprops(nucleiImage,'Centroid'); % Segment the image and return the centre of mass of each region
centroids = cat(1, statsCentroid.Centroid); % Convert result into an arry
figure; imshow(nucleiImage); hold on; % Plot the nuclei image
plot(centroids(:,1), centroids(:,2), 'b*'); hold off; % Overlay the region centres on the image of the nuclei
% Display the number of regions found by region props
title(sprintf('Number of regions: %d',length(centroids)));

