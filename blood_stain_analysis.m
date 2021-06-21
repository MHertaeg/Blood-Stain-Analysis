function [max_location,max_mid_ratio,result] = blood_stain_analysis(image)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Written by Michael Hertaeg

% Image processing algorithm to analyse blood stain patterns for
% diagnostics

% See readme file and publication for details

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% image can be scaled to decrease analysis time for large images.
% Results are not significantly effected by resolution
scale = 0.9;

% Reads in image
A = imread(image);

%Scales image proportionally to "scale" variable
resized = imresize( A , scale ); 

% defines vector of rotations
rot = linspace(0,180,10);

%Loops through each rotation
for u = 1:length(rot)-1
    clearvars -except resized u rot A max_location max_mid_ratio
    
% Rotate image
    if u ~=1
    B = imrotate(resized,rot(u));
    else
        B = resized;
    end

% Select blue channel
c3 = B(:,:,3);

% sets blank pixels outside of image to be white.
[size_vert,size_horiz,kk] = size(B);
for i = 1:size_vert
   for j = 1:size_horiz
        if c3(i,j) == 0
            c3(i,j) = 255;
        end
   end
end


sumrow_hor = zeros(1,size_horiz);
sumrow2 = zeros(1,size_vert);

% Thresholds image
BW = im2bw(c3,0.6);

% Creates vectors summing the horizontal and vertical pixels
for i = 1:size_horiz
    sumrow_hor(i) = sum(BW(:,i));
end
for i = 1:size_vert
    sumrow_vert(i) = sum(BW(i,:));
end

% Find midpoint of stain on horizontal axis
midpoint_horiz = find(sumrow_hor == min(sumrow_hor),1);
%midpoint horiz defined by widest point
condi2 = 0;

% Finds upper and lower boundaries
for i = 1:length(sumrow_vert)
   
    if sumrow_vert(i) < max(sumrow_vert)*0.95 & condi2 <1 & exist('lower2')~=1
        probe(i) = 1;
      
        lower2 = i;
    end
    if sumrow_vert(i) == min(sumrow_vert)
        condi2 = 1;
    end
    
    if condi2 ==1 && sumrow_vert(i) > max(sumrow_vert)*0.95 
        probe(i) = 1;
        condi2 =0;
        upper2 = i;
        break
    end
    
end


range = round(midpoint_horiz/20);

upper2 = min([upper2+range*2,size_vert]);
lower2 = max([lower2-range*2,1]);
cropped = c3(lower2:upper2,midpoint_horiz-range:midpoint_horiz+range);

averow2 = zeros(1,((upper2-(lower2))));
for jjj = 1:((upper2-(lower2)))
    averow2(jjj) = mean(cropped(jjj,:));
end


drift = zeros(1,((upper2)-(lower2)));
for jjj = 1:((upper2)-(lower2))
  
    drift(jjj) = (averow2(1)-averow2(end))/(-length(averow2))*(jjj-1) + averow2(1);
end

% Removes linear drift
Average_inten = -(averow2-drift);


smoothed_intensity = smooth(Average_inten,size_vert/400);
scaled_inten_sm = smoothed_intensity/max(smoothed_intensity);

condi = 0;


for i = 1:length(scaled_inten_sm)
   
    if scaled_inten_sm(i) > 0.4 & condi <1 & exist('lower')~=1
        probe(i) = 1;
       
        
        lower = i;
    end
    if scaled_inten_sm(i) == 1
        condi = 1;
    end
    
    if scaled_inten_sm(i) < 0.4 & condi ==1
        probe(i) = 1;
        condi =0;
        upper = i;
        break
    end
    
end

Maxpoint = find(scaled_inten_sm==1,1);
Midpoint = (upper+lower)/2;

stain_radius = (upper-lower)/2;

max_mid_ratio(u) = mean(scaled_inten_sm(round(Midpoint)-range:round(Midpoint)+range));
max_location(u) = (abs(Maxpoint-Midpoint))/stain_radius;
end
% max_location
max_location = mean(max_location);
% max_mid_ratio
max_mid_ratio = mean(max_mid_ratio);

% Predicts result. These values were predicive in our experience but
% different applications may require adjustments

if max_location >0.8 & max_mid_ratio < 0.8
    result = 'Negative';
elseif max_location < 0.8 & max_mid_ratio > 0.8
     result = 'Positive';
else
    result = 'Error';
end

end
