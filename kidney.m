% Main script to process a single CT-Scan image and visualize each step.
% 1. Image Acquisition
% Load a single image slice. Replace 'sample_image.png' with your file path.
% This example uses a placeholder.
img = imread('/MATLAB Drive/1.3.46.670589.33.1.63700780314615924400001.4656503585389726474.png');
% Convert to grayscale if the image is in color
if size(img, 3) == 3
   img_gray = rgb2gray(img);
else
   img_gray = img;
end
% Display the original image
figure;
imshow(img_gray);
title('1. Original CT-Scan Image');
% 2. Preprocessing
% Apply contrast stretching
img_contrast = imadjust(img_gray);
% Display the image after contrast stretching
figure;
imshow(img_contrast);
title('2. After Contrast Stretching');
% Apply a median filter to reduce noise
img_filtered = medfilt2(img_contrast);
% Display the image after median filtering
figure;
imshow(img_filtered);
title('3. After Median Filtering');
% 3. Segmentation
% Use automatic thresholding to create a binary mask
level = graythresh(img_filtered);
img_binary = imbinarize(img_filtered, level);
% Display the binary image after thresholding
figure;
imshow(img_binary);
title('4. After Thresholding (Binary Mask)');
% Perform morphological erosion to remove small noise or objects
se_erode = strel('disk', 2);
img_eroded = imerode(img_binary, se_erode);
% Display the image after erosion
figure;
imshow(img_eroded);
title('5. After Morphological Erosion');
% Perform morphological dilation to connect parts of the stone and fill holes
se_dilate = strel('disk', 2);
img_segmented = imdilate(img_eroded, se_dilate);
% Display the final segmented mask
figure;
imshow(img_segmented);
title('6. Final Segmented Kidney Stone');
% 4. Volume Estimation (ROI Highlight)
% Overlay the segmented mask on the original image to highlight the stone.
% Convert the original image to a format for overlay
img_rgb = cat(3, img_gray, img_gray, img_gray);
% Create a red mask for the stone area
red_mask = img_segmented;
red_mask(:,:,2:3) = 0; % Set green and blue channels to zero
% Combine the original image and the red mask
highlighted_img = img_rgb;
highlighted_img(repmat(img_segmented, [1, 1, 3])) = 255; % Or a specific color value
% Display the original image with the highlighted stone
figure;
imshow(highlighted_img);
title('7. Original Image with Highlighted ROI');






% Main script to process a single CT-Scan image and segment the kidney stone.

% 1. Image Acquisition
% Load your uploaded image file.
img = imread('1.3.46.670589.33.1.63700700749865510700001.5062181202000819812.png');

% Convert to grayscale if the image is in color
if size(img, 3) == 3
    img_gray = rgb2gray(img);
else
    img_gray = img;
end

% 2. Preprocessing
% Apply a median filter to reduce noise, as described in the paper.
img_filtered = medfilt2(img_gray, [3 3]);

% 3. Segmentation
% Use a high threshold to isolate the brightest regions (stone and bones).
threshold_value = 200;
img_binary = img_filtered > threshold_value;

% Perform morphological operations to refine the mask.
% Remove small objects (noise) that are less than 100 pixels.
img_filtered_size = bwareaopen(img_binary, 100);

% 4. Find the properties of all remaining objects
stats = regionprops(img_filtered_size, 'Area', 'BoundingBox', 'Eccentricity');
all_areas = [stats.Area];

if ~isempty(all_areas)
    % Find objects that are likely kidney stones based on shape and size
    % A stone is typically a round or elliptical object, so its eccentricity
    % should be less than 0.9.
    eccentricities = [stats.Eccentricity];
    
    % Filter out objects that are too elongated (high eccentricity)
    valid_objects_idx = find(eccentricities < 0.9);

    if ~isempty(valid_objects_idx)
        % From the valid objects, find the largest one by area
        valid_stats = stats(valid_objects_idx);
        valid_areas = [valid_stats.Area];
        [~, largest_valid_idx] = max(valid_areas);
        
        % Get the properties of the final selected object
        area_pixels = valid_stats(largest_valid_idx).Area;
        boundingBox = valid_stats(largest_valid_idx).BoundingBox;

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
        disp('No kidney stone detected based on shape and size properties.');
    end
else
    disp('No objects were detected in the image.');
end


% Main script to process a single CT-Scan image and segment the kidney stone.

% 1. Image Acquisition
% Load your uploaded image file.
img = imread('1.3.46.670589.33.1.63700700749865510700001.5062181202000819812.png');

% Convert to grayscale if the image is in color
if size(img, 3) == 3
    img_gray = rgb2gray(img);
else
    img_gray = img;
end

% 2. Preprocessing
% Apply a median filter to reduce noise, as described in the paper.
img_filtered = medfilt2(img_gray, [3 3]);

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


