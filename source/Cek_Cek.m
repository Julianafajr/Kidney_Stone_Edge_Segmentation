clc;
clear;
close all;

% --- 1. Load images (slices) ---
image_files = {
    '1.3.46.670589.33.1.63726694992278881800001.4786973521425870422.png'
    '1.3.46.670589.33.1.63726694992438891000001.5614093351577235839.png'
    
    
};

num_slices = length(image_files);
roi_masks = cell(num_slices, 1);     
gt_masks = cell(num_slices, 1);      
pred_masks = cell(num_slices, 1);    
areas = zeros(num_slices, 1);        

TP_all = zeros(num_slices, 1);
FP_all = zeros(num_slices, 1);
FN_all = zeros(num_slices, 1);

% --- 2. Spatial resolution ---
pixel_spacing = 0.75;           
slice_thickness = 5.0;          
pixel_area = pixel_spacing^2;

% --- 3. Loop through each image ---
for i = 1:num_slices
    I = imread(image_files{i});
    I = mat2gray(I);  

    % --- ROI ginjal manual ---
    figure; imshow(I, []); title(['Draw ROI for kidney in slice ', num2str(i)]);
    roi_mask = roipoly(); close;
    roi_masks{i} = roi_mask;

    % --- GT batu ginjal manual ---
    figure; imshow(I, []); title(['Draw GT for stone in slice ', num2str(i)]);
    gt_mask = roipoly(); close;
    gt_masks{i} = gt_mask;

    % --- Preprocessing ---
    I_eq = adapthisteq(I);
    I_filt = medfilt2(I_eq, [3 3]);
    
    % --- Visualisasi hasil preprocessing ---
    figure;
    subplot(1,2,1); imshow(I, []); title(['Original Slice ', num2str(i)]);
    subplot(1,2,2); imshow(I_filt, []); title(['Preprocessed Slice ', num2str(i)]);

    % Histogram untuk membantu threshold
    figure;
    imhist(I_filt);
    title(['Histogram Slice ', num2str(i)]); %tambahan

    % --- Intensitas rata-rata di dalam ROI ginjal ---
    I_roi = I_filt .* roi_mask;
    vals = I_roi(roi_mask);
    mean_roi = mean(vals);
    std_roi = std(vals);

    % --- Threshold adaptif (lebih sensitif) ---
    %level = mean_roi + 1.0 * std_roi; % awal: 1.5
    %bw = imbinarize(I_filt, level);
    
    % --- Threshold sangat longgar (target: recall tinggi)
    level1 = mean_roi + 0.5 * std_roi; % longgar
    level2 = mean_roi + 0.8 * std_roi; % sedang

    bw1 = imbinarize(I_filt, level1);
    bw2 = imbinarize(I_filt, level2);
    bw_adaptive = imbinarize(I_filt, 'adaptive', 'Sensitivity', 0.45);

    bw = bw1 | bw2 | bw_adaptive;  % kombinasi agresif

    % --- Jika gagal, threshold lebih longgar ---
    if sum(bw(:)) < 10
        level = mean_roi + 0.6 * std_roi;
        bw = imbinarize(I_filt, level);
    end

    % --- Gabung dengan metode lain jika perlu ---
    if sum(bw(:)) < 10
        otsu_level = graythresh(I_filt);
        bw_otsu = imbinarize(I_filt, otsu_level);
        bw_adapt = imbinarize(I_filt, 'adaptive', 'Sensitivity', 0.4);
        bw = bw | bw_otsu | bw_adapt;
    end

    % --- Morfologi ---
    bw = imopen(bw, strel('disk', 1));
    bw = imclose(bw, strel('disk', 3));
    bw = imfill(bw, 'holes');
    bw = imdilate(bw, strel('disk', 2)); % tambahan untuk recall

    % --- Filter berdasarkan ukuran area batu ---
    bw = bwareafilt(bw, [1, 30000]);

    % --- Masking ROI ginjal ---
    pred_mask = bw & roi_mask;
    pred_masks{i} = pred_mask;

    % --- Hitung Area & Evaluasi ---
    areas(i) = sum(pred_mask(:)) * pixel_area;
    TP = sum(pred_mask(:) & gt_mask(:));
    FP = sum(pred_mask(:) & ~gt_mask(:));
    FN = sum(~pred_mask(:) & gt_mask(:));

    TP_all(i) = TP;
    FP_all(i) = FP;
    FN_all(i) = FN;

    % --- Visualisasi segmentasi ---
    figure;
    imshow(I_eq, []); hold on;
    visboundaries(roi_mask, 'Color', 'g', 'LineWidth', 1.5);
    visboundaries(gt_mask, 'Color', 'b', 'LineWidth', 1.5);
    visboundaries(pred_mask, 'Color', 'r', 'LineWidth', 1.5);
    title(['Segmentasi Slice ', num2str(i)]);
    legend({'ROI Ginjal','GT Batu Ginjal','Prediksi'}, 'TextColor','w');
end

% --- 4. Estimasi Volume ---
volumes_mm3 = areas * slice_thickness;
total_volume_mm3 = sum(volumes_mm3);
total_volume_ml = total_volume_mm3 / 1000;

% --- 5. Evaluasi Keseluruhan ---
total_TP = sum(TP_all);
total_FP = sum(FP_all);
total_FN = sum(FN_all);

precision = total_TP / (total_TP + total_FP + eps);
recall    = total_TP / (total_TP + total_FN + eps);
f1_score  = 2 * (precision * recall) / (precision + recall + eps);

% --- 6. Tampilkan Hasil ---
fprintf('\n=== Hasil Estimasi dan Evaluasi ===\n');
fprintf('Total Volume Batu Ginjal: %.2f mm3\n', total_volume_ml);
fprintf('Precision: %.2f%%\n', precision * 100);
fprintf('Recall   : %.2f%%\n', recall * 100);
fprintf('F1 Score : %.2f%%\n', f1_score * 100);

% --- 7. Visualisasi Semua Slice ---
figure;
for i = 1:num_slices
    subplot(1, num_slices, i);
    I = imread(image_files{i});
    imshow(I, []);
    hold on;
    visboundaries(roi_masks{i}, 'Color', 'g');
    visboundaries(gt_masks{i}, 'Color', 'b');
    visboundaries(pred_masks{i}, 'Color', 'r');
    title(['Slice ', num2str(i)]);
end
sgtitle('Ringkasan Segmentasi Seluruh Slice');
