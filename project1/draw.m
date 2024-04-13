% Read an image
I = imread('finaloutput.png');

% Convert the image to grayscale
I_gray = rgb2gray(I);

% Apply Gaussian smoothing
I_smoothed = imgaussfilt(I_gray, 2); % Standard deviation of 2 for smoothing

% Use the Canny edge detector
BW = edge(I_smoothed, 'Canny');

% Close gaps in the edges
BW_closed = bwmorph(BW, 'close');

% Invert the image
BW_inverted = imcomplement(BW_closed);

% Create a new image with black lines on white background
final_image = uint8(255 * BW_inverted);

% Display the final image
imshow(final_image);
title('Final Image');