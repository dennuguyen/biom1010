rgbImage = imread('ExampleFingerprint.jpg');
grayImage = rgb2gray(rgbImage); % Convert to grayscale
figure; imshow(grayImage); % Display Image

croppedImage = imcrop(grayImage); % Interactive cropping tool

se1 = strel('disk',5); % Structuring element
filteredImage = imtophat(croppedImage,se1); % Top-hat filtering
figure; imshow(filteredImage); % Display figure
figure; imhist(filteredImage); % Display histogram

contrastImage = imadjust(filteredImage); % Contrast stretching
figure; imshow(contrastImage); % Display figure
figure; imhist(contrastImage); % Display histogram

otsuLevel = graythresh(contrastImage); % Get Otsu threshold 
bwImage = im2bw(contrastImage,otsuLevel); % Convert to black and white
figure; imshow(bwImage); % Display black and white image

se2 = strel('disk',5); % Structuring element
temp = imclose(bwImage,se2); % Morphological closing
filledImage = imopen(temp,se2); % Morphological opening
figure; imshow(filledImage); % Display image

smoothedImage = medfilt2(filledImage,[3,3]); % Median filtering
figure; imshow(smoothedImage); % Display image

% Remove residual areas having less than 50 pixels
contactImage = bwareaopen(smoothedImage,50);
figure; imshow(contactImage); % Display image

statsArea = regionprops(contactImage,'Area');
area = struct2array(statsArea);

