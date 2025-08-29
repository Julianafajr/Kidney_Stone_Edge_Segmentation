% Main script to process a single CT-Scan image and segment the kidney stone.
% 1. Image Acquisition
% Load your uploaded image file.
img = imread('/MATLAB Drive/1.3.46.670589.33.1.63700700750188529100001.5659992199131706929.png');
% Convert to grayscale if the image is in color
if size(img, 3) == 3
   img_gray = rgb2gray(img);
else
   img_gray = img;
end
%display image
figure;
imshow(img_gray);
title('1.Original CT-Scan Image');
%display the histogram of the orignal image
figure;
imhist(img_gray);
title('Histogram of Image')
% 2. Preprocessing
% Apply a median filter to reduce noise.
img_filtered = medfilt2(img_gray, [3 3]);
%display image
figure;
imshow(img_filtered);
title('2. Image after median filter');

% 3. Segmentation
% Use a high threshold to isolate the brightest regions (stone and bones).
threshold_value = 200;
img_binary = img_filtered > threshold_value;
% Perform morphological operations to refine the mask.
% Close the small holes within the segmented region.
se_close = strel('disk', 5);
img_segmented = imclose(img_binary, se_close);
% Remove any small noise objects (less than 50 pixels).
img_segmented = bwareaopen(img_segmented, 50);
%display segmented binary image
figure;
imshow(img_segmented);
title('3.Segmented binary mask');

% 4. Find the properties of all remaining objects
stats = regionprops(img_segmented, 'Area', 'BoundingBox', 'Eccentricity');
all_areas = [stats.Area];
if ~isempty(all_areas)
   % Filter out objects that are too large (likely the spine or other major bone)
   max_area_threshold = 5000; % You may need to adjust this value
   size_filtered_stats = stats([stats.Area] < max_area_threshold);
   if ~isempty(size_filtered_stats)
       % Now, find the most circular object from the remaining list.
       % This is a more robust approach than a fixed eccentricity value.
       eccentricities = [size_filtered_stats.Eccentricity];
      
       % Sort objects by eccentricity in ascending order (closest to 0 is a perfect circle)
       [sorted_ecc, sorted_idx] = sort(eccentricities, 'ascend');
      
       % The most circular object is the first one in the sorted list
       final_idx = sorted_idx(1);
      
       % Get the properties of the final selected object
       area_pixels = size_filtered_stats(final_idx).Area;
       boundingBox = size_filtered_stats(final_idx).BoundingBox;

       % 5. Volume Estimation
       % The paper assumes a slice thickness of 5 mm.
       slice_thickness = 5; % in mm
       estimated_volume = area_pixels * slice_thickness;
      
       fprintf('Estimated volume of the kidney stone: %.2f mm^3\n', estimated_volume);
       % Display the original image with the red bounding box
       figure;
       imshow(img_gray);
       title('Kidney Stone Highlighted with Red Bounding Box');
       hold on;
      
       % Draw the bounding box in red ('r').
       rectangle('Position', boundingBox, 'EdgeColor', 'r', 'LineWidth', 2);
       hold off;
   else
       disp('All objects were too large to be a kidney stone.');
   end
else
   disp('No objects were detected in the image.');
end


