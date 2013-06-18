clear all;
close all;
fileFolder = 'J:/LaptopSource/Matlab/EE243/HW1/Data/Day3_Scene28_Camera22/';
 
dirOutput=dir(fullfile(fileFolder, '*.bmp'));
%sort in descending order of modified date
% the key frames are extracted first
% the background frames are extracted later
[Y fileSortIndex] = sort({dirOutput.date});   
fileNames={dirOutput(fileSortIndex).name}';

bgImageNum = 199;
% key frames with moving people
image_names = char(fileNames( end-bgImageNum-199:2:end-bgImageNum ));

% Processing one image and only save the final results for report
selectFrames = 1:1:size(image_names,1);

aviobj = avifile('example.avi','compression','None');
fig=figure;
for i=selectFrames
    im_RGB = double(imread( fullfile(fileFolder, image_names(i,:) ) ));
    imshow(im_RGB./256);
    
    F = getframe(fig);
    aviobj = addframe( aviobj,F );
end
close(fig);
aviobj = close(aviobj);