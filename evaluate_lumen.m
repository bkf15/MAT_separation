labels = get_labels();
concordant_samples = zeros(1, 2);
for i = 1:size(labels.p1_labels, 1)
	if labels.p1_labels(i) == labels.p2_labels(i) && labels.p1_labels(i) == labels.p3_labels(i)
		concordant_samples = [concordant_samples; labels.p1_labels(i) labels.slide_names(i)];
	end
end
concordant_samples = concordant_samples(2:end, :);
%ignore 'other' samples
fea_samples = concordant_samples(concordant_samples(:, 1) == 'Flat Epithelial', :);
columnar_samples = concordant_samples(concordant_samples(:, 1) == 'Columnar', :);
adh_samples = concordant_samples(concordant_samples(:, 1) == 'ADH', :);
normal_samples = concordant_samples(concordant_samples(:, 1) == 'Normal Duct', :);

for i = 1:size(normal_samples,1)
	im_path = char(normal_samples(i, 2));
	im_path = char(pwd + "/all_slides/" + im_path(1:strfind(im_path,'_')-1) + "/" + im_path);
	figure; imshow(imread(im_path));
end
%select the best examples of lumen segmentation
normal_samples = normal_samples([6 14 29 56 62], :);
fea_samples = fea_samples([1 3 4 5 6], :);
columnar_samples = columnar_samples([3 6 10 13 14], :);
adh_samples = adh_samples([1 3 5 6 8], :);

normal_scores = [];
fea_scores = [];
columnar_scores = [];
adh_scores = [];

for i = 1:size(normal_samples, 1)
	disp(string(i));
	im_path = char(normal_samples(i, 2));
	%im_path = char(pwd + "/all_slides/" + im_path(1:strfind(im_path,'_')-1) + "/" + im_path);
	[~, lumen_scores] = compute_statistics(im_path, 1);
	normal_scores = [normal_scores ; lumen_scores];
	
	im_path = char(fea_samples(i, 2));
	%im_path = char(pwd + "/all_slides/" + im_path(1:strfind(im_path,'_')-1) + "/" + im_path);
	[~, lumen_scores] = compute_statistics(im_path, 1);
	fea_scores = [fea_scores ; lumen_scores];
	
	im_path = char(columnar_samples(i, 2));
	%im_path = char(pwd + "/all_slides/" + im_path(1:strfind(im_path,'_')-1) + "/" + im_path);
	[~, lumen_scores] = compute_statistics(im_path, 1);
	columnar_scores = [columnar_scores ; lumen_scores];
	
	im_path = char(adh_samples(i, 2));
	%im_path = char(pwd + "/all_slides/" + im_path(1:strfind(im_path,'_')-1) + "/" + im_path);
	[~, lumen_scores] = compute_statistics(im_path, 1);
	adh_scores = [adh_scores ; lumen_scores];
end

%euclidian distance as the variance metric
%ribbon symmetry
variance_matrix = zeros(4);

variance_matrix(1, 1) = norm(normal_scores(:, 1) - (ones(size(normal_scores(:,1)))*mean(normal_scores(:,1))));
variance_matrix(2, 2) = norm(adh_scores(:, 1) - (ones(size(adh_scores(:,1)))*mean(adh_scores(:,1))));
variance_matrix(3, 3) = norm(fea_scores(:, 1) - (ones(size(fea_scores(:,1)))*mean(fea_scores(:,1))));
variance_matrix(4, 4) = norm(columnar_scores(:, 1) - (ones(size(columnar_scores(:,1)))*mean(columnar_scores(:,1))));

variance_matrix(1, 2) = norm(normal_scores(:, 1) - adh_scores(:, 1));
variance_matrix(1, 3) = norm(normal_scores(:, 1) - fea_scores(:, 1));
variance_matrix(1, 4) = norm(normal_scores(:, 1) - columnar_scores(:, 1));
variance_matrix(2, 3) = norm(adh_scores(:, 1) - fea_scores(:, 1));
variance_matrix(2, 4) = norm(adh_scores(:, 1) - columnar_scores(:, 1));
variance_matrix(3, 4) = norm(fea_scores(:, 1) - columnar_scores(:, 1));

%plot the variance matrix
t = num2cell(variance_matrix);
t = cellfun(@num2str, t, 'UniformOutput', false); 
X = [1 2 3 4 2 3 4 3 4 4];
Y = [1 1 1 1 2 2 2 3 3 4];
t = [t(1,1) t(1,2) t(1,3) t(1,4) t(2,2) t(2,3) t(2,4) t(3,3) t(3,4) t(4,4)]; 
figure; imagesc(variance_matrix); text(X, Y, t, 'HorizontalAlignment', 'Center'); colorbar;
title(sprintf('Lumen ribbon within/between class variance\n Normal               ADH               FEA                 Columnar'));
ylabel('Columnar             FEA             ADH               Normal');

%now taper symmetry
%euclidian distance
variance_matrix = zeros(4);

variance_matrix(1, 1) = norm(normal_scores(:, 2) - (ones(size(normal_scores(:,2)))*mean(normal_scores(:,2))));
variance_matrix(2, 2) = norm(adh_scores(:, 2) - (ones(size(adh_scores(:,2)))*mean(adh_scores(:,2))));
variance_matrix(3, 3) = norm(fea_scores(:, 2) - (ones(size(fea_scores(:,2)))*mean(fea_scores(:,2))));
variance_matrix(4, 4) = norm(columnar_scores(:, 2) - (ones(size(columnar_scores(:,2)))*mean(columnar_scores(:,2))));

variance_matrix(1, 2) = norm(normal_scores(:, 2) - adh_scores(:, 2));
variance_matrix(1, 3) = norm(normal_scores(:, 2) - fea_scores(:, 2));
variance_matrix(1, 4) = norm(normal_scores(:, 2) - columnar_scores(:, 2));
variance_matrix(2, 3) = norm(adh_scores(:, 2) - fea_scores(:, 2));
variance_matrix(2, 4) = norm(adh_scores(:, 2) - columnar_scores(:, 2));
variance_matrix(3, 4) = norm(fea_scores(:, 2) - columnar_scores(:, 2));

%plot the variance matrix
t = num2cell(variance_matrix);
t = cellfun(@num2str, t, 'UniformOutput', false); 
X = [1 2 3 4 2 3 4 3 4 4];
Y = [1 1 1 1 2 2 2 3 3 4];
t = [t(1,1) t(1,2) t(1,3) t(1,4) t(2,2) t(2,3) t(2,4) t(3,3) t(3,4) t(4,4)]; 
figure; imagesc(variance_matrix); text(X, Y, t, 'HorizontalAlignment', 'Center'); colorbar;
title(sprintf('Lumen taper within/between class variance\n Normal               ADH               FEA                 Columnar'));
ylabel('Columnar             FEA             ADH               Normal');


%separation
%euclidian distance
variance_matrix = zeros(4);

variance_matrix(1, 1) = norm(normal_scores(:, 3) - (ones(size(normal_scores(:,3)))*mean(normal_scores(:,3))));
variance_matrix(2, 2) = norm(adh_scores(:, 3) - (ones(size(adh_scores(:,3)))*mean(adh_scores(:,3))));
variance_matrix(3, 3) = norm(fea_scores(:, 3) - (ones(size(fea_scores(:,3)))*mean(fea_scores(:,3))));
variance_matrix(4, 4) = norm(columnar_scores(:, 3) - (ones(size(columnar_scores(:,3)))*mean(columnar_scores(:,3))));

variance_matrix(1, 2) = norm(normal_scores(:, 3) - adh_scores(:, 3));
variance_matrix(1, 3) = norm(normal_scores(:, 3) - fea_scores(:, 3));
variance_matrix(1, 4) = norm(normal_scores(:, 3) - columnar_scores(:, 3));
variance_matrix(2, 3) = norm(adh_scores(:, 3) - fea_scores(:, 3));
variance_matrix(2, 4) = norm(adh_scores(:, 3) - columnar_scores(:, 3));
variance_matrix(3, 4) = norm(fea_scores(:, 3) - columnar_scores(:, 3));

%plot the variance matrix
t = num2cell(variance_matrix);
t = cellfun(@num2str, t, 'UniformOutput', false); 
X = [1 2 3 4 2 3 4 3 4 4];
Y = [1 1 1 1 2 2 2 3 3 4];
t = [t(1,1) t(1,2) t(1,3) t(1,4) t(2,2) t(2,3) t(2,4) t(3,3) t(3,4) t(4,4)]; 
figure; imagesc(variance_matrix); text(X, Y, t, 'HorizontalAlignment', 'Center'); colorbar;
title(sprintf('Lumen separation within/between class variance\n Normal               ADH               FEA                 Columnar'));
ylabel('Columnar             FEA             ADH               Normal');
