%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% #STAR - CORNER DETECTION %%
%% Shayan Fazeli - 91102171 %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%preparing the script:
clear all;
close all;
clc;

%reading the image:
disp('please enter the number of star you want to use as input.');
disp('for example, if you want star1.bmp, type in 1 and press enter');
num=input('My selection is ');
mystar = imread(['star', sprintf('%d', num),'.bmp']);

%%%%%%%%%%%%%%%%%%%%%%%%%
%PERFORMING HARRIS CORNER DETECTION
%well, what we do here is not "exactly" harris corner detection, but
%it's a variant and the ideas are exactly the same and gotten by that
%article.


%okay, first we need to have a gaussian type filter
%but it's derivative, in order to obtain the Ix and Iy
%needed for harris corner detector, we need to define it.
%the problem is there is no function equivalent for this purpose
%which was predefined in matlab, so we start by defining the
%gaussian derivative filter.

%window_size is going to determine the height and width of gaussian
%window we are willing to convolve our image with. it will be scaled
%automatically with the image size (number of rows to be exact).
%the value of the parameter has been set as a function, which it's
%relation itself has been deducted empirically and it might be 
%a good idea to change it for a specific scale that has not been
%taken into consideration in this script.
window_size = ceil(5*(size(mystar,1)/400));

%now, we need to define a meshgrid, based on which we define our gaussians:
[X, Y] = meshgrid(-(floor(window_size/2)):floor(window_size/2),...
        -(floor(window_size/2)):floor(window_size/2));
    
%standard deviation will be used in our gaussians. by trial and error
%it's been deduced that the value 1 works out to be the best choice.
standard_deviation = 1;
%what is the derivative of a gaussian? we use it here, to create horizonal
%and vertical gradient filters:
horizontal_gradient_filter = X .* exp(-(X.^2+Y.^2)/(2*standard_deviation^2));
vertical_gradient_filter = Y .* exp(-(X.^2+Y.^2)/(2*standard_deviation^2));
%and finally, the gaussian filter itself:
gaussian_filter = exp(-(X.^2+Y.^2)/(2*standard_deviation^2));

%now, we convolve these gaussian filters, in other words we are
%using a "weighted" gradient filter:
Ix = conv2(horizontal_gradient_filter, double(rgb2gray(mystar)));
Iy = conv2(vertical_gradient_filter, double(rgb2gray(mystar)));

%now, smoothening up the things:
Ix2 = conv2(gaussian_filter, Ix.^2);
Iy2 = conv2(gaussian_filter, Iy.^2);
Ixy = conv2(gaussian_filter, Ix.*Iy);

%seeing harris corner detection algorithm's details in
%wikipedia, it says that the cost function defined in the article
%sometimes gets tough to work with. so, we use this cost function
%which is more common. first, we define a small positive integer, called
%the epsilon:
epsilon = 0.00001;

%now for "each" pixel, we define an energy value:
M = 2*(Ix2.*Iy2-(Ixy.^2))./(Ix2+Iy2+epsilon);
%taking a copy of it:
M2 = M;

%in this stage, we show the energy value by which we are going to decide 
%where are the potential corner points and their whereabouts:
imagesc(M); title('energy function');

disp('script is paused, press enter to continue..');
pause;

%after pause, it's time to deal with the energy map.
%first, we need to put a threshold in the field, and throw
%the pixels with energy less than that threshold out as garbage.

%empirically, it's been concluded that a threshold, scalable
%with the size and about 30 percent of the maximum value of M
%(based on imhist and estimating the variance and investigating
%the distribution of the values), is going to work good in this case.
threshold = 0.3 * (size(mystar,1)/400) * max(M(:));
%compute the indices of pixels with energy less than that:
indices = M < threshold;
%throwing them out:
M(indices)=0;



%now, i try to find maximums, then burn them down...
corners = [];
cutter_half_window = floor(10); %this could be scaled, yet 
%for the range we're working on i didnt think it was necessary, though
%it is possible that in a specific case this be better changed.
%25 times finding the maximums:
for i = 1:25
    [xmax1, ymax1] = find(M == (max(M(:))));
    xmax = xmax1(1);
    ymax = ymax1(1);
    %in above lines we have found the first maximum we could have found
    %now, concatenating corners with it:
    corners = [corners; [xmax, ymax]];
    %wiping it's whereabouts clean so we don't get stuck in there:
    M(max(1,xmax-cutter_half_window):min(size(M,1),xmax+cutter_half_window),...
        max(1,ymax-cutter_half_window):min(size(M,2),ymax+cutter_half_window))=0;
end



%now, just in case, purging the corners array:
%first, we define new corners:
new_corners = [];
%first, add the first corner to the array:
new_corners = [new_corners; corners(1,:)];
%now for the other corners...
i = 2;
%this flag will be used in the loop:
flag = 0;
%purge threshold means the minimum distance
%two corner points are allowed to have. i haven't worked
%long on it's scaling mechanism, but if the scale being used
%be abnormal, it could simply be changed and the code still works.
%for the scales close to our examples, this simple two mode threshold
%will do.
purge_threshold = 20;
if (size(mystar,1)<100)
    purge_threshold=10;
end

%so, sweeping the corners array:
while i <= size(corners,1)
    %check whether or not there is a point in the array that
    %is taken into consideration before as another corner point:
    for j=1:size(new_corners,1)
        if (i>=0 && i<=size(corners,1))
            if(norm(new_corners(j,:)-corners(i,:))<purge_threshold)
                %if so, do this:
                 new_corners(j,:) = 0.5*(new_corners(j,:)+corners(i,:));
                i = i + 1;
                flag = 1;
            end
        else
            break;
        end
    end
    if flag == 1
        flag = 0;
    else
        %if not, add it up, continue:
        new_corners = [new_corners; corners(i,:)];
        i = i + 1;
    end
end
%saving the changes:
corners = new_corners;

%now, it's time to show the star and it's corners:
imshow(mystar); title('the star and its corners'); hold on;

temp_corners = [];
%showing them using red Xs..
for i =1:size(corners,1)
    if(~((corners(i,1)<10 && corners(i,2)<10) || ...
            ((corners(i,1)>size(M,1)-10) && (corners(i,2) > size(M,2)-10)) ||...
            ((corners(i,1)>size(M,1)-10) && corners(i,2)<10) || ...
            ((corners(i,1)<10 && (corners(i,2) > size(M,2)-10)))))
        plot(corners(i,2),corners(i,1),'rx');
        temp_corners=[temp_corners; corners(i,:)];
    end
end


%now, trying to know them.
%first, we need to compute the center:
corners = temp_corners;
center = mean(corners);

%now we have to compute the distance of the points, from the center:
differences = corners - repmat(center,[size(corners,1),1]);
distances = zeros(size(differences,1),1);
for i = 1:size(differences,1)
    distances(i,1) = norm(differences(i,:));
end

%now, we have two clusters, one, the inner circle, the points closer
%to the center:
cluster1=zeros(10,1);
[useless, index] = min(distances);
distances(index,1)=Inf;
cluster1(index,1)=1;
[useless, index] = min(distances);
distances(index,1)=Inf;
cluster1(index,1)=1;
[useless, index] = min(distances);
distances(index,1)=Inf;
cluster1(index,1)=1;
[useless, index] = min(distances);
distances(index,1)=Inf;
cluster1(index,1)=1;
[useless, index] = min(distances);
distances(index,1)=Inf;
cluster1(index,1)=1;

%now the cluster2:
cluster2 = ~cluster1;

pcluster1 = corners(find(cluster1),:);
pcluster2 = corners(find(cluster2),:);

%now we get one point, compute the distance between this point
%that was chosen from cluster1, and all the points from
%the cluster2, that will give us the edge of the star. now,
%using that size, we can manage to connect them.
[cpoints1, cpoints2] = find_where_to_go(pcluster1, pcluster2);


%now that we have the points, and what they are connected to, let us 
%draw:
shapeInserter = vision.ShapeInserter('Shape', 'Lines', 'BorderColor', 'Custom', 'CustomBorderColor', [0,0,255], 'LineWidth', 1);

for i = 1:5
   point1 = pcluster1(i,:);
   point2 = cpoints1(i,:);
   point3 = cpoints2(i,:);
   mystar = step(shapeInserter, mystar, int32(cat(2,flip(point1), flip(point2))));
   mystar = step(shapeInserter, mystar, int32(cat(2,flip(point1), flip(point3))));
end

disp('Script is paused, continue to see the final result.')
pause;
close all;
figure;
title('FINAL');
imshow(mystar);
hold on;
for i = 1:size(corners,1)
   plot(corners(i,2), corners(i,1), 'rx'); 
end











