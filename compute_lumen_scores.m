function [lumen_score] = compute_lumen_scores(lumen_boundary_scores)
	num_lumen = 3;
	lumen_CC = bwconncomp(lumen_boundary_scores);
	lumen_scores = [];
	for i = 1:lumen_CC.NumObjects
		avg_score = sum(lumen_boundary_scores(lumen_CC.PixelIdxList{i})) / size(lumen_CC.PixelIdxList{i}, 1);
		lumen_scores = [lumen_scores; avg_score i];
	end
	lumen_scores = sortrows(lumen_scores, 'descend');
	if size(lumen_scores, 1) == 0
		lumen_score = 0;
		disp('Zero lumen error');
		return;
	end
	if size(lumen_scores, 1) < num_lumen
		num_lumen = size(lumen_scores, 1);
	end
	lumen_score = sum(lumen_scores(1:num_lumen,1)) / num_lumen;
end