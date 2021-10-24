image = imread('FakeTan.jpg'); % Read RGB image
R = image(:,:,1); % Red component of image
G = image(:,:,2); % Green component of image
B = image(:,:,3); % Blue component of image

condition1 = R > 95 & G > 40 & B > 20; % Condition 1
condition2 = (max(image,[],3)-min(image,[],3)) > 15; % Condition 2
condition3 = abs(R-G) > 15 & R > G & R > B; % Condition 3

% Is skin if all 3 conditions satisfied
isSkinImage = condition1 & condition2 & condition3;
% Plot original image and logical image of skin
figure; 
subplot(2,1,1); imshow(image);
subplot(2,1,2); imshow(isSkinImage);
