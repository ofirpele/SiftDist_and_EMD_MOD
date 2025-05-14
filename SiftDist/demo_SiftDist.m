% SiftDist_demo Demonstrate SiftDist.  
%   This demo computes the sift descriptor on a pair of well known
%   test images. Then it matches them using SiftDist.

clc;

fprintf(1,'The script calls functions from Andrea Vedaldi`s sift code.\n');
fprintf(1,'You need to:\n');
fprintf(1,' 1. download his code from:\n');
  
fprintf(1,'     http://www.vlfeat.org/~vedaldi/assets/sift/binaries/  \n');
fprintf(1,'       ( it contains already compiled binaries, for the source code see: )\n');
fprintf(1,'       ( http://www.vlfeat.org/~vedaldi/assets/sift/versions/ )\n');
fprintf(1,' 2. Unzip the file.\n');
fprintf(1,' 3. Change the path in line 17 to the directory in which Andrea Vedaldi`s sift code is.\n');

%addpath path_to_Andrea_Vedaldi_code;
addpath sift/

% Increase value to decrease the number of matches 
% (and possibly better matches).
% Run time is expected to increase/decrease 
% when DistRatio increase/decrease.
% Note that the returned ratios for inds~=0 
% are always DistRatio.
% If set to -1 all ratios are computed and the running time is full.
% Look at DistRatio param of SiftRatioMatch
% for more details.
DistRatio= 1.25;
NumOrientBins= 16;

I1=imreadbw('img1.ppm');
I2=imreadbw('img3.ppm');

I1=I1-min(I1(:));
I1=I1/max(I1(:));
I2=I2-min(I2(:));
I2=I2/max(I2(:));

fprintf('Computing frames and descriptors (~0.5 minutes).\n');
[frames1,descr1]= sift(I1, 'NumOrientBins', NumOrientBins);
[frames2,descr2]= sift(I2, 'NumOrientBins', NumOrientBins);

fprintf('Matches using SiftDist (~1 minutes with default params).\n');
% Taking the sqrt of descr improves results.
% I still have no theoretical reason for why it improves results.
[inds ratios] = SiftRatioMatch(sqrt(descr1), frames1, sqrt(descr2), frames2, DistRatio);


% Transform to Andrea's Vedaldi's matches format
matches= zeros(2,sum(inds~=0));
i= 1;
for r=1:size(descr1,2)
  if (inds(r)~=0)
    matches(1,i)= r;
    matches(2,i)= inds(r);
    i=i+1;
  end
end

figure; clf;
title('Click on a keypoint to see its match'); axis off;
plotmatches(I1, I2, frames1(1:2,:), frames2(1:2,:), matches, 'Interactive', 2);
drawnow;


