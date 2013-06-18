clear all
close all

SAVE_IMAGE = 1;

data = xlsread( 'Point coordinates' );
data = data';

% calibIndex = 1:size(data,2);
calibIndex = 1:20;
points = data(:,calibIndex);
% load input
% the input is a matrix named 'points'. Each column is a data vector for a point.
% The entries of each data vector is
% point ID - 3D wolrd coordinate - 2D image coordinate

fileFolder = '../Data/Day3_Scene28_Camera22/';
 
dirOutput=dir(fullfile(fileFolder, '*.bmp'));
%sort in descending order of modified date
% the key frames are extracted first
% the background frames are extracted later
[Y fileSortIndex] = sort({dirOutput.date});   
fileNames={dirOutput(fileSortIndex).name}';

% background images: there is 20 frames extracted as background images
bgImageNum = 20;
background_names = char(fileNames(end-bgImageNum+1:end));
% read an image to initialize the background image
background_image = double(imread( fullfile(fileFolder, background_names(end,:) ) ));
figure
imshow(background_image./256)
hold on
% plot the dot representing a corner feature
plot( points(5,:), points(6,:), 'r.', 'MarkerSize', 18 )
plot( points(5,1), points(6,1), 'g.', 'MarkerSize', 18 )
hold off

% Save the scene with feature points selected
if (SAVE_IMAGE)
    % Saving the plot as .fig
    hgsave('Report/plots/select_points')
    % Saving the plot as .esp (color) with 300 dpi resolution
    print -depsc2 -r300 'Report/plots/select_points'
    % Saving the plot as .jpg image
    saveas(gcf,'Report/plots/select_points.jpg')
end

% Camera calibration
ox = 320;
oy = 240;

points(5,:) = points(5,:) - ox;
points(6,:) = points(6,:) - oy;
points(2:4,:) = points(2:4,:).*0.01; % from cm to m

totalFeatNum = size(points,2);
A_mat = zeros( totalFeatNum, 8 );
for i=1:totalFeatNum
    A_mat(i,1) = points(5,i)*points(2,i);
    A_mat(i,2) = points(5,i)*points(3,i);
    A_mat(i,3) = points(5,i)*points(4,i);
    A_mat(i,4) = points(5,i);
    A_mat(i,5) = -points(6,i)*points(2,i);
    A_mat(i,6) = -points(6,i)*points(3,i);
    A_mat(i,7) = -points(6,i)*points(4,i);
    A_mat(i,8) = -points(6,i);
end

[U S V] = svd(A_mat);
% the solution is the eigenvector corresponding to the smallest eigenvalue
v = V(:,end);

r2 = v(1:3);
% find the scaling factor
gamma = 1/sqrt(r2'*r2);
r1 = v(5:7);
alpha = sqrt(r1'*r1)/sqrt(r2'*r2);

r2 = gamma*r2;
r1 = gamma*r1/alpha;
r3 = cross(r1,r2);
t1 = v(8)/alpha*gamma;
t2 = v(4)*gamma;

selectIdx = 1;
if ( points(5,selectIdx)*(r1'*points(2:4,selectIdx) + t1) > 0 )
    r1 = -r1;
    t1 = -t1;
end

if ( points(6,selectIdx)*(r2'*points(2:4,selectIdx) + t2) > 0 )
    r2 = -r2;
    t2 = -t2;
end

R = [r1'; r2'; r3'];
[U S V] = svd(R);
R = U*V';

A_mat = zeros(totalFeatNum,2);
b_vec = zeros(totalFeatNum,1);
for i=1:totalFeatNum
    A_mat(i,1) = points(5,i);
    A_mat(i,2) = points(2:4,i)'*r1 + t1;
    b_vec(i) = -points(5,i)*(points(2:4,i)'*r3);
end
sol = pinv(A_mat)*b_vec;
t3 = sol(1);

% Output
aspect_ratio = alpha

rotation_matrix = R

translation_vector = [t1 t2 t3]'

effective_focal_length = sol(2)

% Reprojection
data(2:4,:) = data(2:4,:).*0.01; % from cm to m
projection = zeros(2,size(data,2));
for i=1:size(data,2)
    c_p = rotation_matrix*data(2:4,i) + translation_vector;
    projection(1,i) = -effective_focal_length*c_p(1)/c_p(3) + ox;
    projection(2,i) = -effective_focal_length/aspect_ratio*c_p(2)/c_p(3) + oy;
end

cooordinate_axes = [0 1 0 0 0 0;0 0 0 1 0 0;0 0 0 0 0 1];
coord_projection = zeros(2,size(cooordinate_axes,2));
for i=1:size(cooordinate_axes,2)
    c_p = rotation_matrix*cooordinate_axes(:,i) + translation_vector;
    coord_projection(1,i) = -effective_focal_length*c_p(1)/c_p(3) + ox;
    coord_projection(2,i) = -effective_focal_length/aspect_ratio*c_p(2)/c_p(3) + oy;
end

figure
imshow(background_image./256)
hold on
% plot the dot representing a corner feature
plot( data(5,:), data(6,:), 'r.', 'MarkerSize', 18 )
plot( projection(1,:), projection(2,:), 'b.', 'MarkerSize', 18 )
plot( coord_projection(1,:), coord_projection(2,:), 'k-','LineWidth',2 )
hold off

% Save the scene with feature points selected
if (SAVE_IMAGE)
    % Saving the plot as .fig
    hgsave('Report/plots/reprojection_points')
    % Saving the plot as .esp (color) with 300 dpi resolution
    print -depsc2 -r300 'Report/plots/reprojection_points'
    % Saving the plot as .jpg image
    saveas(gcf,'Report/plots/reprojection_points.jpg')
end
