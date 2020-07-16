function [scores] = find_high_symmetry_nuclei(boundary_symmetry, slide_im, labels, impath)
	CC = bwconncomp(boundary_symmetry);
	%scores will be Nx2 where N = number of CC's. first column contains
	%that CC's normalized sym score, second column is its index in
	%CC.PixelIdxList
	scores = [];
	for i = 1:CC.NumObjects
		[r, c] = ind2sub(size(boundary_symmetry), CC.PixelIdxList{i});
		%if connected component < 12 pixels, don't count it
		if size(r, 1) < 12
			continue
		end
		%sum symmetry scores on contour and divide by number of pixels on
		%the contour
		sym_total = 0;
		for j = 1:(size(CC.PixelIdxList{i}))
			sym_total = sym_total + boundary_symmetry(r(j), c(j));
		end
		sym_total = sym_total / size(CC.PixelIdxList{i}, 1);
		scores = [sym_total, i; scores];
	end
	%sort the scores so the highest are at the top
	if isempty(scores)
		disp("itsaempty");
	end
	scores = sortrows(scores, 'descend');
		
	%{
	%find the index of the slide in labels, to access doctor labels
	slide_file_name = regexp(impath, '/adh03-2a_');
	slide_index = find(labels.slide_names == impath(slide_file_name+1:end));
	%loop through the sorted scores array and display the nuclei in
	%descending order when a button is pressed 
	f = figure; 
	f.Name = string(impath(slide_file_name+1:end));
	for i = 1:size(scores, 1)
		[r, c] = ind2sub(size(boundary_symmetry), CC.PixelIdxList{scores(i, 2)});
		try
			%draw a red circle around the region of interest
			subplot(1, 2, 1); imagesc(boundary_symmetry(r(1)-50:r(1)+50, c(1)-50:c(1)+50)); colorbar;
			title('Ribbon symmetry');
			slide_im_circ = insertShape(slide_im, 'circle', [c(1), r(1), 40], 'Color', 'r', 'LineWidth', 1);
			subplot(1, 2, 2); imshow(slide_im_circ(r(1)-100:r(1)+100, c(1)-100:c(1)+100, :));
			label_string = [labels.participant_1 " " labels.p1_labels(slide_index) newline labels.participant_2 " " labels.p2_labels(slide_index) newline labels.participant_3 " " labels.p3_labels(slide_index)];
			annotation('textbox', [.01 .5 .3 .3], 'String', string(label_string), 'FitBoxToText', 'on');
			w = waitforbuttonpress;
			clf;
		catch e
			%if an error is thrown, its because the region is too close to
			%the edge to display
			clf;
			%"%Index in position 1 is invalid. Array indices must be
			%positive integers or logical values."
			if strlength(e.message) > 65
				%index = which dimension is out of bounds 
				index = str2num(e.message(19));
				if index == 1
					if c(1) >= 100
						subplot(1, 2, 1); imagesc(boundary_symmetry(1:r(1)+50, c(1)-50:c(1)+50)); colorbar;
						title('Ribbon symmetry');
						slide_im_circ = insertShape(slide_im, 'circle', [c(1), r(1), 40], 'Color', 'r', 'LineWidth', 1);
						subplot(1, 2, 2); imshow(slide_im_circ(1:r(1)+100, c(1)-100:c(1)+100, :));
						label_string = [labels.participant_1 " " labels.p1_labels(slide_index) newline labels.participant_2 " " labels.p2_labels(slide_index) newline labels.participant_3 " " labels.p3_labels(slide_index)];
						annotation('textbox', [.01 .5 .3 .3], 'String', string(label_string), 'FitBoxToText', 'on');
						w = waitforbuttonpress;
					else
						subplot(1, 2, 1); imagesc(boundary_symmetry(1:r(1)+50, 1:c(1)+50)); colorbar;
						title('Ribbon symmetry');
						slide_im_circ = insertShape(slide_im, 'circle', [c(1), r(1), 40], 'Color', 'r', 'LineWidth', 1);
						subplot(1, 2, 2); imshow(slide_im_circ(1:r(1)+100, 1:c(1)+100, :));
						label_string = [labels.participant_1 " " labels.p1_labels(slide_index) newline labels.participant_2 " " labels.p2_labels(slide_index) newline labels.participant_3 " " labels.p3_labels(slide_index)];
						annotation('textbox', [.01 .5 .3 .3], 'String', string(label_string), 'FitBoxToText', 'on');
						w = waitforbuttonpress;
					end
				else
					subplot(1, 2, 1); imagesc(boundary_symmetry(r(1)-50:r(1)+50, 1:c(1)+50)); colorbar;
					title('Ribbon symmetry');
					slide_im_circ = insertShape(slide_im, 'circle', [c(1), r(1), 40], 'Color', 'r', 'LineWidth', 1);
					subplot(1, 2, 2); imshow(slide_im_circ(r(1)-100:r(1)+100, 1:c(1)+100, :));
					label_string = [labels.participant_1 " " labels.p1_labels(slide_index) newline labels.participant_2 " " labels.p2_labels(slide_index) newline labels.participant_3 " " labels.p3_labels(slide_index)];
					annotation('textbox', [.01 .5 .3 .3], 'String', string(label_string), 'FitBoxToText', 'on');
					w = waitforbuttonpress;
				end
			%Index in position 2 exceeds array bounds (must not exceed 1025).'
			else
				index = str2num(e.message(19));
				if index == 1
					subplot(1, 2, 1); imagesc(boundary_symmetry(r(1)-50:end, c(1)-50:c(1)+50)); colorbar;
					title('Ribbon symmetry');
					slide_im_circ = insertShape(slide_im, 'circle', [c(1), r(1), 40], 'Color', 'r', 'LineWidth', 1);
					subplot(1, 2, 2); imshow(slide_im_circ(r(1)-100:end, c(1)-100:c(1)+100, :));
					label_string = [labels.participant_1 " " labels.p1_labels(slide_index) newline labels.participant_2 " " labels.p2_labels(slide_index) newline labels.participant_3 " " labels.p3_labels(slide_index)];
					annotation('textbox', [.01 .5 .3 .3], 'String', string(label_string), 'FitBoxToText', 'on');
					w = waitforbuttonpress;
				else
					subplot(1, 2, 1); imagesc(boundary_symmetry(r(1)-50:r(1)+50, c(1)-50:end)); colorbar;
					title('Ribbon symmetry');
					slide_im_circ = insertShape(slide_im, 'circle', [c(1), r(1), 40], 'Color', 'r', 'LineWidth', 1);
					subplot(1, 2, 2); imshow(slide_im_circ(r(1)-100:r(1)+100, c(1)-100:end, :));
					label_string = [labels.participant_1 " " labels.p1_labels(slide_index) newline labels.participant_2 " " labels.p2_labels(slide_index) newline labels.participant_3 " " labels.p3_labels(slide_index)];
					annotation('textbox', [.01 .5 .3 .3], 'String', string(label_string), 'FitBoxToText', 'on');
					w = waitforbuttonpress;
				end
			end
			clf;
		end
	end
	%}
end