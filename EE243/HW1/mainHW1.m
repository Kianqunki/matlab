clear all;
close all;
fileFolder = 'J:/LaptopSource/Matlab/EE243/HW1/Data/Day3_Scene28_Camera22/';
 
dirOutput=dir(fullfile(fileFolder, '*.bmp'));
%sort in descending order of modified date
% the key frames are extracted first
% the background frames are extracted later
[Y fileSortIndex] = sort({dirOutput.date});   
fileNames={dirOutput(fileSortIndex).name}';

% background images: there is 199 frames extracted as background images
bgImageNum = 20;
background_names = char(fileNames(end-bgImageNum+1:end));
% % check the images
% for i=1:20:size(background_names,1)
%     im = imread( fullfile(fileFolder, background_names(i,:) ) );
%     figure
%     imshow(im)
% end

% read a dummy image to initialize the background image
background_image = double(imread( fullfile(fileFolder, background_names(1,:) ) ));
for i=1:bgImageNum
    im = double(imread( fullfile(fileFolder, background_names(i,:) ) ));
    background_image = background_image + im;
end
background_image_RGB = background_image/bgImageNum;
[rowNum, colNum, colorDepth] = size(background_image);
background_image = rgb2hsv( background_image_RGB );
% % check background image
% figure, imshow(background_image./256);

bgImageNum = 199;
% key frames with moving people
image_names = char(fileNames( end-bgImageNum-199:2:end-bgImageNum ));
% % check the images
% for i=1:10:size(image_names,1)
%     im = imread( fullfile(fileFolder, image_names(i,:) ) );
%     figure
%     imshow(im)
% end

% the threshold to determine that the color is not from background
thres = 30;
thresH = 0.35;
thresS = 0.25;
thresV = 70;

% human component size
humanSizeThres = 3000;

% Image window code
% 101: foreground image
% 102: current frame
% 103: background image
% 104: "cleaned" foreground image
% 105: "cleaned" foreground image 2
% 106: labeled foreground image
% 107: labeled foreground image

figure(102), imshow(background_image_RGB./256);
% Saving the plot as .fig
hgsave('Report/plots/background')
% Saving the plot as .esp (color) with 300 dpi resolution
print -depsc2 -r300 'Report/plots/background'
% Saving the plot as .jpg image
saveas(gcf,'Report/plots/background.jpg')

% Processing one image and save all intermediate steps fpr report
selectFrames=4;

for i=selectFrames
    im_RGB = double(imread( fullfile(fileFolder, image_names(i,:) ) ));
    % convert to HSV color space
    im = rgb2hsv( im_RGB );
    
    % check current image
    figure(103), imshow(im_RGB./256);
    % Saving the plot as .fig
    hgsave('Report/plots/input_image')
    % Saving the plot as .esp (color) with 300 dpi resolution
    print -depsc2 -r300 'Report/plots/input_image'
    % Saving the plot as .jpg image
    saveas(gcf,'Report/plots/input_image.jpg')
    
    % perform background subtraction in RGB color space
    foreground_image_RGB = zeros(rowNum,colNum);
    foreground_image_RGB = ( abs(im_RGB(:,:,1) - background_image_RGB(:,:,1)) > thres ... 
        | abs(im_RGB(:,:,2) - background_image_RGB(:,:,2)) > thres ...
        | abs(im_RGB(:,:,3) - background_image_RGB(:,:,3)) > thres );
    
    figure(101), imshow(foreground_image_RGB);
    % Saving the plot as .fig
    hgsave('Report/plots/foreground_RGB')
    % Saving the plot as .esp (color) with 300 dpi resolution
    print -depsc2 -r300 'Report/plots/foreground_RGB'
    % Saving the plot as .jpg image
    saveas(gcf,'Report/plots/foreground_RGB.jpg')
    
    % by delating first, we connect and isolate the larger components
    clean_foreground = foreground_image_RGB;
    clean_foreground = bwmorph( clean_foreground, 'dilate', 2 );
    clean_foreground = bwmorph( clean_foreground, 'erode', 2 );
    % find the number of connected components in foreground
    [label_foreground, compNum] = bwlabeln( clean_foreground, 8 );
    % remove the small connected components
    for k=1:compNum
        [r,c] = find( label_foreground == k );
        if ( length(r) < humanSizeThres )
            clean_foreground(r,c) = 0;
        end
    end
    clean_foreground_RGB = logical(clean_foreground);
    
    figure(104), imshow(clean_foreground_RGB);
    % Saving the plot as .fig
    hgsave('Report/plots/clean_foreground_RGB')
    % Saving the plot as .esp (color) with 300 dpi resolution
    print -depsc2 -r300 'Report/plots/clean_foreground_RGB'
    % Saving the plot as .jpg image
    saveas(gcf,'Report/plots/clean_foreground_RGB.jpg')
    
    % RGB is generally good but having a problem:
    % false positive against shadows. Therefore, using HSV to imrpove on it
    foreground_image = zeros(rowNum,colNum);
    foreground_image = ( (abs(im(:,:,1) - background_image(:,:,1)) > thresH ... 
        | abs(im(:,:,2) - background_image(:,:,2)) > thresS ...
        | abs(im(:,:,3) - background_image(:,:,3)) > thresV) );
    
    figure(101), imshow(foreground_image);
    % Saving the plot as .fig
    hgsave('Report/plots/foreground_RGB')
    % Saving the plot as .esp (color) with 300 dpi resolution
    print -depsc2 -r300 'Report/plots/foreground_HSV'
    % Saving the plot as .jpg image
    saveas(gcf,'Report/plots/foreground_HSV.jpg')
    
    % perform erode and dilate operations to remove "noises"
    % by delating first, we connect and isolate the larger components
    clean_foreground = foreground_image;
    clean_foreground = bwmorph( clean_foreground, 'dilate', 2 );
    clean_foreground = bwmorph( clean_foreground, 'erode', 2 );
    % find the number of connected components in foreground
    [label_foreground, compNum] = bwlabeln( clean_foreground, 8 );
    % remove the small connected components
    for k=1:compNum
        [r,c] = find( label_foreground == k );
        if ( length(r) < humanSizeThres )
            clean_foreground(r,c) = 0;
        end
    end
    clean_foreground_HSV = logical(clean_foreground);
    
    figure(105), imshow(clean_foreground_HSV);
    % Saving the plot as .fig
    hgsave('Report/plots/clean_foreground_HSV')
    % Saving the plot as .esp (color) with 300 dpi resolution
    print -depsc2 -r300 'Report/plots/clean_foreground_HSV'
    % Saving the plot as .jpg image
    saveas(gcf,'Report/plots/clean_foreground_HSV.jpg')    
    
    % Use the cleaned foreground as a mask
    clean_foreground2 = clean_foreground_RGB & clean_foreground_HSV;
    clean_foreground2 = bwmorph( clean_foreground2, 'dilate', 5 );
    clean_foreground2 = bwmorph( clean_foreground2, 'erode', 5 );    
    
    figure(104), imshow(clean_foreground_RGB & clean_foreground_HSV);
    figure(105), imshow(clean_foreground2);
    
    % find the number of connected components in foreground
    [label_foreground2, compNum] = bwlabeln( clean_foreground2, 8 );
    
    % check the labels
    figure(106), imshow( label_foreground2./compNum);
    % Saving the plot as .fig
    hgsave('Report/plots/segmentation')
    % Saving the plot as .esp (color) with 300 dpi resolution
    print -depsc2 -r300 'Report/plots/segmentation'
    % Saving the plot as .jpg image
    saveas(gcf,'Report/plots/segmentation.jpg') 
    
    temp_im = im_RGB./256;
    temp_im2 = im_RGB./256;
    for k=1:compNum
        [r, c] = find( label_foreground2 == k  );
        temp_im(r,c,k) = 1.0;
        for l=1:length(r)
            temp_im2(r(l),c(l),k) = 1.0;
        end
%         keyboard;
    end
    
    % check the labels
    figure(108), imshow( temp_im );
    % Saving the plot as .fig
    hgsave('Report/plots/final')
    % Saving the plot as .esp (color) with 300 dpi resolution
    print -depsc2 -r300 'Report/plots/final'
    % Saving the plot as .jpg image
    saveas(gcf,'Report/plots/final', 'jpg')
    
    % check the labels
    figure(109), imshow( temp_im2 );
    % Saving the plot as .fig
    hgsave('Report/plots/final2')
    % Saving the plot as .esp (color) with 300 dpi resolution
    print -depsc2 -r300 'Report/plots/final2'
    % Saving the plot as .jpg image
    saveas(gcf,'Report/plots/final2', 'jpg')

end

% Processing one image and only save the final results for report
selectFrames = 5:5:size(image_names,1);

for i=selectFrames
    im_RGB = double(imread( fullfile(fileFolder, image_names(i,:) ) ));
    % convert to HSV color space
    im = rgb2hsv( im_RGB );
    
    % perform background subtraction in RGB color space
    foreground_image_RGB = zeros(rowNum,colNum);
    foreground_image_RGB = ( abs(im_RGB(:,:,1) - background_image_RGB(:,:,1)) > thres ... 
        | abs(im_RGB(:,:,2) - background_image_RGB(:,:,2)) > thres ...
        | abs(im_RGB(:,:,3) - background_image_RGB(:,:,3)) > thres );
    % by delating first, we connect and isolate the larger components
    clean_foreground = foreground_image_RGB;
    clean_foreground = bwmorph( clean_foreground, 'dilate', 2 );
    clean_foreground = bwmorph( clean_foreground, 'erode', 2 );
    % find the number of connected components in foreground
    [label_foreground, compNum] = bwlabeln( clean_foreground, 8 );
    % remove the small connected components
    for k=1:compNum
        [r,c] = find( label_foreground == k );
        if ( length(r) < humanSizeThres )
            clean_foreground(r,c) = 0;
        end
    end
    clean_foreground_RGB = logical(clean_foreground);
    
    % RGB is generally good but having a problem:
    % false positive against shadows. Therefore, using HSV to imrpove on it
    foreground_image = zeros(rowNum,colNum);
    foreground_image = ( (abs(im(:,:,1) - background_image(:,:,1)) > thresH ... 
        | abs(im(:,:,2) - background_image(:,:,2)) > thresS ...
        | abs(im(:,:,3) - background_image(:,:,3)) > thresV) );
    % perform erode and dilate operations to remove "noises"
    % by delating first, we connect and isolate the larger components
    clean_foreground = foreground_image;
    clean_foreground = bwmorph( clean_foreground, 'dilate', 2 );
    clean_foreground = bwmorph( clean_foreground, 'erode', 2 );
    % find the number of connected components in foreground
    [label_foreground, compNum] = bwlabeln( clean_foreground, 8 );
    % remove the small connected components
    for k=1:compNum
        [r,c] = find( label_foreground == k );
        if ( length(r) < humanSizeThres )
            clean_foreground(r,c) = 0;
        end
    end
    clean_foreground_HSV = logical(clean_foreground);
    
%     fore1 = zeros(rowNum,colNum);
%     fore2 = zeros(rowNum,colNum);
%     fore3 = zeros(rowNum,colNum);
%     fore1 = logical(clean_foreground) & abs(im(:,:,1) - background_image(:,:,1)) > thresH;
%     fore2 = logical(clean_foreground) & abs(im(:,:,2) - background_image(:,:,2)) > thresS;
%     fore3 = logical(clean_foreground) & abs(im(:,:,3) - background_image(:,:,3)) > thresV;
%     foreground_image = (fore1 | fore2) | fore3;
%     figure(111),imshow(fore1);
%     figure(112),imshow(fore2);
%     figure(113),imshow(fore3);            
    
    % Use the cleaned foregrounds as a mask
    clean_foreground2 = clean_foreground_RGB & clean_foreground_HSV;
    clean_foreground2 = bwmorph( clean_foreground2, 'dilate', 5 );
    clean_foreground2 = bwmorph( clean_foreground2, 'erode', 5 );    
    
%     figure(104), imshow(clean_foreground_RGB & clean_foreground_HSV);
%     figure(105), imshow(clean_foreground2);
    
    % find the number of connected components in foreground
    [label_foreground2, compNum] = bwlabeln( clean_foreground2, 8 );
    
    temp_im = im_RGB./256;
    temp_im2 = im_RGB./256;
    for k=1:compNum
        [r, c] = find( label_foreground2 == k  );
        temp_im(r,c,k) = 1.0;
        for l=1:length(r)
            temp_im2(r(l),c(l),k) = 1.0;
        end
    end
    
    filename1 = sprintf( 'Report/plots/final_frame%d', i);
    filename2 = sprintf( 'Report/plots/final2_frame%d', i);

    % check the labels
    figure(108), imshow( temp_im );
    % Saving the plot as .fig
    hgsave(filename1)
    % Saving the plot as .esp (color) with 300 dpi resolution
    print( '-depsc2', '-r300', filename1 );
    % Saving the plot as .jpg image
    saveas(gcf,filename1, 'jpg')

    % check the labels
    figure(109), imshow( temp_im2 );
    % Saving the plot as .fig
    hgsave(filename2)
    % Saving the plot as .esp (color) with 300 dpi resolution
    print( '-depsc2', '-r300', filename2 );
    % Saving the plot as .jpg image
    saveas(gcf,filename2, 'jpg')
    
    i
end