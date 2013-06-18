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

% the threshold to determine that the color is not from background
thres = 30;
thresH = 0.35;
thresS = 0.25;
thresV = 70;
interThresUpp = 0.6;
interThresLow = 0.1;

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

% Kalman filter initial parameters
R=[ 25 5;
    5 25];
H= [1 0 0 0;
    0 1 0 0];
Q= 10*eye(4);
% std of initial position error is 10 pixels
P= 100*eye(4);
dt=2; % because we process every two frames
F= [1 0 dt 0;
    0 1 0 dt;
    0 0 1 0;
    0 0 0 1];
kalmanInit = 0;
% use each Kalman filter to independently to track one target
% Maximum 5 target
% Each filter has x: state vector, P: covariance matrix, template: the last
% detected target
maxFilterNum = 5;
filters = repmat(struct('x', [], 'P', 100*eye(4), 'template', zeros(rowNum, colNum), 'frame', zeros(1,4) ), maxFilterNum, 1);
% descriptors store which filter is initialized, and tracking which object
descriptors = [];

% Processing one image and only save the final results for report
selectFrames = 1:2:size(image_names,1);
savePos = zeros(2,size(image_names,1),maxFilterNum);

for i=selectFrames
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Feature detection step
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
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
        
    % Use the cleaned foregrounds as a mask
    clean_foreground2 = clean_foreground_RGB & clean_foreground_HSV;
    figure(104), imshow(clean_foreground2);
    
    clean_foreground2 = bwmorph( clean_foreground2, 'dilate', 5 );
    clean_foreground2 = bwmorph( clean_foreground2, 'erode', 5 );    
        
    % find the number of connected components in foreground
    [label_foreground2, compNum] = bwlabeln( clean_foreground2, 8 );
    % remove the small connected components
    for k=1:compNum
        [r,c] = find( label_foreground2 == k );
        if ( length(r) < humanSizeThres )
            clean_foreground2(r,c) = 0;
        end
    end
    [label_foreground2, compNum] = bwlabeln( clean_foreground2, 8 );
    
%     temp_im = im_RGB./256;
% %     temp_im2 = im_RGB./256;
%     for k=1:compNum
%         [r, c] = find( label_foreground2 == k  );
%         temp_im(min(r):max(r),min(c):max(c),k) = 1.0;
% %         for l=1:length(r)
% %             temp_im2(r(l),c(l),k) = 1.0;
% %         end
%     end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Initialize the first filters
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    if ( isempty(descriptors) )
        for k=1:compNum
            descriptors = [descriptors k];
            
            [r, c] = find( label_foreground2 == k  );
            posR = (min(r)+max(r))/2;
            posC = (min(c)+max(c))/2;
            filters(k).x = [posR posC 0 0]';
            for l=1:length(r)
                filters(k).template(r(l),c(l)) = 1.0;
            end
            filters(k).frame = [min(r) max(r) min(c) max(c)];
            
%             % check
%             figure(109), imshow(filters(k).template)
        end
        
        savePos(:,i,k) = filters(k).x(1:2);
        % skip the rest
        continue;
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Filter propagation and data association
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    % flag for data association
    % 0: for not matched, -1: for ignore, ID: for matching
    dataAssocFlag = zeros(1,compNum);
    z = repmat(struct('descriptors', [], 'meas', zeros(2, 1), 'template', zeros(rowNum, colNum)), compNum, 1);
    for j=1:length(descriptors)
        prevState = filters(j).x;

        % propagate the state
        filters(j).x = F*prevState;
        diffPos = floor(filters(j).x(1:2) - prevState(1:2));

        % propagate the covariance
        filters(j).P = F*filters(j).P*F' + Q;

        % propagate the template
        minRow = max( filters(j).frame(1)+diffPos(1), 0 );
        maxRow = min( filters(j).frame(2)+diffPos(1), rowNum );
        minCol = max( filters(j).frame(3)+diffPos(2), 0 );
        maxCol = min( filters(j).frame(4)+diffPos(2), colNum );
        propTemplate = zeros(rowNum,colNum);
        propTemplate(minRow:maxRow,minCol:maxCol) = ...
            filters(j).template(filters(j).frame(1):filters(j).frame(2),filters(j).frame(3):filters(j).frame(4));
        currTemplateSize = length(find(propTemplate>0));

        for k=find( dataAssocFlag == 0 )
            currObj = zeros( rowNum, colNum );
            [r, c] = find( label_foreground2 == k  );
            currObjSize = length(r);
            for l=1:length(r)
                currObj(r(l),c(l)) = 1.0;
            end

            interObjTemplate = logical(currObj) & logical(propTemplate);
            interObjTemplateSize = length(find(interObjTemplate>0));
            interObjTemplateSize/currTemplateSize
            interObjTemplateSize/currObjSize

            % only one detected object for each component
            figure(108), imshow(filters(j).template(filters(j).frame(1):filters(j).frame(2),filters(j).frame(3):filters(j).frame(4)))
            figure(109), imshow(propTemplate)
            figure(110), imshow(currObj)
            figure(111), imshow(interObjTemplate)

            if ( interObjTemplateSize/currTemplateSize > interThresUpp && interObjTemplateSize/currObjSize > interThresUpp )
%             if ( interObjTemplateSize/currTemplateSize > interThresUpp )
                % if the overlapping is significant for both templates, found new measurement
                dataAssocFlag(k) = j;

                z(k).descriptors = j;
                z(k).template = currObj;
                z(k).frame = [min(r) max(r) min(c) max(c)];
                z(k).meas = [(min(r) + max(r))/2 (min(c) + max(c))/2 ]';

                break;
            elseif ( interObjTemplateSize/currTemplateSize < interThresLow && interObjTemplateSize/currObjSize < interThresLow )
                % if not, others' measurement
            elseif ( interObjTemplateSize/currTemplateSize > interThresUpp && interObjTemplateSize/currObjSize < interThresLow )
                % two people close together and detected as one component
            else
                % if somewhere in the middle, bad segmentation, ignore it.
%                 dataAssocFlag(k) = -1;
            end % end if
        end % end for k in all components
    end % end for j in all Kalman filters

    % if there are components not matched, flag = 0, initialize new
    % filter
    idx = find( dataAssocFlag == 0 );
    for k=idx
        newId = length(descriptors)+1;
        descriptors = [descriptors newId];

        [r, c] = find( label_foreground2 == k );
        posR = (min(r)+max(r))/2;
        posC = (min(c)+max(c))/2;
        filters(newId).x = [posR posC 0 -20]';
        for l=1:length(r)
            filters(newId).template(r(l),c(l)) = 1.0;
        end
        filters(newId).frame = [min(r) max(r) min(c) max(c)];

        %             % check
        %             figure(109), imshow(filters(k).template)
        savePos(:,i,newId) = filters(newId).x(1:2);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Filter update
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    for k=find( dataAssocFlag > 0 )
        if ( ~isempty(z(k).descriptors) )
            % the ID of the Kalman filter
            j = z(k).descriptors;
            filters(j).x
            % predicted measurement
            z_k1 = H*filters(j).x;
            % measurement innovation
            r_k1 = z(k).meas - z_k1;
            S = H*filters(j).P*H' + R;
            % Kalman gain
            K = filters(j).P*H'*inv(S);
            % update
            filters(j).x = filters(j).x + K*r_k1;
            filters(j).P = filters(j).P - K*H*filters(j).P;
            filters(j).template = z(k).template;
            filters(j).frame = z(k).frame;
            filters(j).x
            
%             pause;
        end
    end
    
    % save the tracked positions
    for j=1:length(descriptors)
        savePos(:,i,j) = filters(j).x(1:2);
    end
    
    i
end

im_RGB = double(imread( fullfile(fileFolder, image_names(selectFrames(end),:) ) ));
figure(200), imshow(im_RGB./256);
hold on
plot( savePos(2,selectFrames,1), savePos(1,selectFrames,1), 'g*' );
plot( savePos(2,selectFrames,2), savePos(1,selectFrames,2), 'b*' );
% plot( savePos(2,selectFrames,3), savePos(1,selectFrames,3), 'r*' );
hold off
% filename2 = 'tracking_result2';
% % Saving the plot as .fig
% hgsave(filename2)
% % Saving the plot as .esp (color) with 300 dpi resolution
% print( '-depsc2', '-r300', filename2 );
% % Saving the plot as .jpg image
% saveas(gcf,filename2, 'jpg')