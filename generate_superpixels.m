%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [BW] = generate_superpixels(impath, sk_save_path)
	RGB_Im = imread(impath);
	%disp(f(f_idx).name);
	% Create Superpixels and find stats %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	[L,n] = superpixels(uint8(RGB_Im),300);
	BW = boundarymask(L);
	%imwrite(BW, sk_save_path);
	imshow(imoverlay(RGB_Im, BW, 'cyan'));

%imshow(imoverlay(RGB_Im, BW, 'cyan'));
%CC = regionprops(L);
%XY = cat(1,CC.Centroid);
%TRI = delaunay(XY(:,1), XY(:,2));
%triplot(TRI, XY(:,1), XY(:,2));

%%%%% Get histogram information for all superpixels %%%%%
% Select superpixels with some rbg constraints %%%%%%%%%%%
%idx = label2idx(L);
%numRows = size(RGB_Im,1);
%numCols = size(RGB_Im,2);

%BW = zeros(1025, 1025);
%for k = 1:N
%    redIdx = idx{k};
%    greenIdx = idx{k} + numRows*numCols;
%    blueIdx = idx{k} + 2*numRows*numCols;
    
    %%26-bin histogram 
%    [yRed,x] = imhist(RGB_Im(redIdx), 26);
%    [yBlue,x] = imhist(RGB_Im(blueIdx), 26);
%    [yGreen,x] = imhist(RGB_Im(greenIdx), 26);

%    Sp_hist(k).sp_id = k;
%    Sp_hist(k).red = yRed;
%    Sp_hist(k).blue = yBlue;
%    Sp_hist(k).green = yGreen; 
%end
%imshow(imoverlay(RGB_Im, BW, 'green'));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
