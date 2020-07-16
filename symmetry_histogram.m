%create a histogram of size 2k to store the count of different
%	symmetry scores on the symmetry image 
function [hist] = symmetry_histogram(sym_image, orig_boundaries, k)
	hist = zeros(1, 2*k);
	%symmetry scores are between 0 and 1, with steps of 1/2k
	%propagate an array of the possible symmetry scores to compare to those
	%present in the image
	sym_values = zeros(1, 2*k);
	v = 0;
	for i = 2:2*k
		v = v + (1 / (2*k));
		sym_values(i) = v;
	end
	
	%loop over symmetry image, if sym_image(i, j) = sym_values(k), then 
	% hist(k) = hist(k) + 1
	
	%get the connected components of the contour image 
	%	NOTE**** this computation is done ~3 times throughout the entire
	%	run of the algorithm. This slows it down and will be fixed in the
	%	future
	CC = bwconncomp(orig_boundaries);
	%loop over each connected component
	for i = 1:CC.NumObjects
		%loop over the pixels in the current CC
		for j = 1:(size(CC.PixelIdxList{i}))
			[y, x] = ind2sub(size(sym_image), CC.PixelIdxList{i}(j));
			sym = sym_image(y, x);
			%note: this for loop exists because find(sym_scores == sym) was
			% not working for certain values. probably a rounding error...
			for k = 1:size(sym_values, 2)
				if ((sym >= (sym_values(k) - 0.01)) && (sym <= (sym_values(k) + 0.01)))
					hist(k) = hist(k) + 1;
				end
			end
		end
	end
	%normalize the histogram
	hist = hist / sum(hist);
end
