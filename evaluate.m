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
%smallest set is fea, with 6 concordant examples

%each slide image has 6 scores, col 1 = nuc_ribbon, col 2 = nuc_taper, col 3 = nuc_separation
%	recall that these scores are averages over the top num_nuclei symmetryic
%	nuclei for each measurement. IE, fea_scores(1, 1) is the average ribbon
%	symmetry score over the top num_nuclei ribbon symmetric nuclei. 
%	the last 3 columns are the lumen ribbon, taper, and separation scores.
%	these are the scores for only the top 1 most symmetric lumen region
fea_scores = [];
columnar_scores = [];
adh_scores = [];
normal_scores = [];

%how many of the top scored nuclei will be considered in the evaluation
num_nuclei = 100;
for i = 1:min([size(fea_samples, 1), size(normal_samples, 1), size(adh_samples, 1), size(columnar_samples, 1)])
	disp(i + " fea");
	[t_nuc,t_lumen] = compute_statistics(fea_samples(i, 2), num_nuclei);
	fea_scores = [fea_scores ; t_nuc t_lumen];
	
	disp(i + " col");
	[t_nuc, t_lumen] = compute_statistics(columnar_samples(i, 2), num_nuclei);
	columnar_scores = [columnar_scores ; t_nuc t_lumen];
	
	disp(i + " adh");
	[t_nuc, t_lumen] = compute_statistics(adh_samples(i, 2), num_nuclei);
	adh_scores = [adh_scores ; t_nuc t_lumen];
	
	disp(i + " normal");
	[t_nuc, t_lumen] = compute_statistics(normal_samples(i, 2), num_nuclei);
	normal_scores = [normal_scores ; t_nuc t_lumen];
end

%create the norm/variance matrix, matrix is KxK, where K = number of
%classes (4, norm/adh/fea/columnar). Diagonal is the within class standard
%deviation for a score. other points are distances from one classes score
%to another
%order is normal, adh, fea, columnarl
%first, ribbon symmetry
%euclidian distance as the variance metric
variance_matrix = zeros(4);

variance_matrix(1, 1) = norm(normal_scores(:, 1) - (ones(size(normal_scores(:,1)))*mean(normal_scores(:,1))));
variance_matrix(2, 2) = norm(columnar_scores(:, 1) - (ones(size(columnar_scores(:,1)))*mean(columnar_scores(:,1))));
variance_matrix(3, 3) = norm(fea_scores(:, 1) - (ones(size(fea_scores(:,1)))*mean(fea_scores(:,1))));
variance_matrix(4, 4) = norm(adh_scores(:, 1) - (ones(size(adh_scores(:,1)))*mean(adh_scores(:,1))));

variance_matrix(1, 2) = norm(normal_scores(:, 1) - columnar_scores(:, 1));
variance_matrix(1, 3) = norm(normal_scores(:, 1) - fea_scores(:, 1));
variance_matrix(1, 4) = norm(normal_scores(:, 1) - adh_scores(:, 1));
variance_matrix(2, 3) = norm(columnar_scores(:, 1) - fea_scores(:, 1));
variance_matrix(2, 4) = norm(columnar_scores(:, 1) - adh_scores(:, 1));
variance_matrix(3, 4) = norm(fea_scores(:, 1) - adh_scores(:, 1));


%angle between the vectors as a variance metric
%{
variance_matrix(1, 1) = acosd(dot(normal_scores(:,1),ones(size(normal_scores(:,1)))*mean(normal_scores(:,1)))/(norm(normal_scores(:,1))*norm(ones(size(normal_scores(:,1)))*mean(normal_scores(:,1)))));
variance_matrix(2, 2) = acosd(dot(adh_scores(:,1),ones(size(adh_scores(:,1)))*mean(adh_scores(:,1)))/(norm(adh_scores(:,1))*norm(ones(size(adh_scores(:,1)))*mean(adh_scores(:,1)))));
variance_matrix(3, 3) = acosd(dot(fea_scores(:,1),ones(size(fea_scores(:,1)))*mean(fea_scores(:,1)))/(norm(fea_scores(:,1))*norm(ones(size(fea_scores(:,1)))*mean(fea_scores(:,1)))));
variance_matrix(1, 1) = acosd(dot(columnar_scores(:,1),ones(size(columnar_scores(:,1)))*mean(columnar_scores(:,1)))/(norm(columnar_scores(:,1))*norm(ones(size(columnar_scores(:,1)))*mean(columnar_scores(:,1)))));

variance_matrix(1,2) = acosd(dot(normal_scores(:,1),adh_scores(:,1)) / (norm(normal_scores(:,1)) * norm(adh_scores(:,1))));
variance_matrix(1,3) = acosd(dot(normal_scores(:,1),fea_scores(:,1)) / (norm(normal_scores(:,1)) * norm(fea_scores(:,1))));
variance_matrix(1,4) = acosd(dot(normal_scores(:,1),columnar_scores(:,1)) / (norm(normal_scores(:,1)) * norm(columnar_scores(:,1))));
variance_matrix(2,3) = acosd(dot(adh_scores(:,1),fea_scores(:,1)) / (norm(adh_scores(:,1)) * norm(fea_scores(:,1))));
variance_matrix(2,4) = acosd(dot(adh_scores(:,1),columnar_scores(:,1)) / (norm(adh_scores(:,1)) * norm(columnar_scores(:,1))));
variance_matrix(3,4) = acosd(dot(fea_scores(:,1),columnar_scores(:,1)) / (norm(fea_scores(:,1)) * norm(columnar_scores(:,1))));
%}

%hausdorff distance
%{
variance_matrix(1,1) = HausdorffDist(normal_scores(:,1),ones(size(normal_scores(:,1)))*mean(normal_scores(:,1)));
variance_matrix(2,2) = HausdorffDist(adh_scores(:,1),ones(size(adh_scores(:,1)))*mean(adh_scores(:,1)));
variance_matrix(3,3) = HausdorffDist(fea_scores(:,1),ones(size(fea_scores(:,1)))*mean(fea_scores(:,1)));
variance_matrix(4,4) = HausdorffDist(columnar_scores(:,1),ones(size(columnar_scores(:,1)))*mean(columnar_scores(:,1)));

variance_matrix(1,2) = HausdorffDist(normal_scores(:,1), adh_scores(:,1));
variance_matrix(1,3) = HausdorffDist(normal_scores(:,1), fea_scores(:,1));
variance_matrix(1,4) = HausdorffDist(normal_scores(:,1), columnar_scores(:,1));
variance_matrix(2,3) = HausdorffDist(adh_scores(:,1), fea_scores(:,1));
variance_matrix(2,4) = HausdorffDist(adh_scores(:,1), columnar_scores(:,1));
variance_matrix(3,3) = HausdorffDist(fea_scores(:,1), columnar_scores(:,1));
%}
t = num2cell(variance_matrix);
t = cellfun(@num2str, t, 'UniformOutput', false); 
X = [1 2 3 4 2 3 4 3 4 4];
Y = [1 1 1 1 2 2 2 3 3 4];
t = [t(1,1) t(1,2) t(1,3) t(1,4) t(2,2) t(2,3) t(2,4) t(3,3) t(3,4) t(4,4)]; 
figure; imagesc(variance_matrix); text(X, Y, t, 'HorizontalAlignment', 'Center'); colorbar;
title(sprintf('Ribbon Symmetry Class Variances\n Normal               Columnar               FEA                 ADH'));
ylabel('ADH             FEA             Columnar               Normal');

%now taper symmetry
%euclidian distance
variance_matrix = zeros(4);

variance_matrix(1, 1) = norm(normal_scores(:, 2) - (ones(size(normal_scores(:,2)))*mean(normal_scores(:,2))));
variance_matrix(2, 2) = norm(columnar_scores(:, 2) - (ones(size(columnar_scores(:,2)))*mean(columnar_scores(:,2))));
variance_matrix(3, 3) = norm(fea_scores(:, 2) - (ones(size(fea_scores(:,2)))*mean(fea_scores(:,2))));
variance_matrix(4, 4) = norm(adh_scores(:, 2) - (ones(size(adh_scores(:,2)))*mean(adh_scores(:,2))));

variance_matrix(1, 2) = norm(normal_scores(:, 2) - columnar_scores(:, 2));
variance_matrix(1, 3) = norm(normal_scores(:, 2) - fea_scores(:, 2));
variance_matrix(1, 4) = norm(normal_scores(:, 2) - adh_scores(:, 2));
variance_matrix(2, 3) = norm(columnar_scores(:, 2) - fea_scores(:, 2));
variance_matrix(2, 4) = norm(columnar_scores(:, 2) - adh_scores(:, 2));
variance_matrix(3, 4) = norm(fea_scores(:, 2) - adh_scores(:, 2));

%angle between the vectors as a variance metric
%{
variance_matrix(1, 1) = acosd(dot(normal_scores(:,2),ones(size(normal_scores(:,2)))*mean(normal_scores(:,2)))/(norm(normal_scores(:,2))*norm(ones(size(normal_scores(:,2)))*mean(normal_scores(:,2)))));
variance_matrix(2, 2) = acosd(dot(adh_scores(:,2),ones(size(adh_scores(:,2)))*mean(adh_scores(:,2)))/(norm(adh_scores(:,2))*norm(ones(size(adh_scores(:,2)))*mean(adh_scores(:,2)))));
variance_matrix(3, 3) = acosd(dot(fea_scores(:,2),ones(size(fea_scores(:,2)))*mean(fea_scores(:,2)))/(norm(fea_scores(:,2))*norm(ones(size(fea_scores(:,2)))*mean(fea_scores(:,2)))));
variance_matrix(1, 1) = acosd(dot(columnar_scores(:,2),ones(size(columnar_scores(:,2)))*mean(columnar_scores(:,2)))/(norm(columnar_scores(:,2))*norm(ones(size(columnar_scores(:,2)))*mean(columnar_scores(:,2)))));

variance_matrix(1,2) = acosd(dot(normal_scores(:,2),adh_scores(:,2)) / (norm(normal_scores(:,2)) * norm(adh_scores(:,2))));
variance_matrix(1,3) = acosd(dot(normal_scores(:,2),fea_scores(:,2)) / (norm(normal_scores(:,2)) * norm(fea_scores(:,2))));
variance_matrix(1,4) = acosd(dot(normal_scores(:,2),columnar_scores(:,2)) / (norm(normal_scores(:,2)) * norm(columnar_scores(:,2))));
variance_matrix(2,3) = acosd(dot(adh_scores(:,2),fea_scores(:,2)) / (norm(adh_scores(:,2)) * norm(fea_scores(:,2))));
variance_matrix(2,4) = acosd(dot(adh_scores(:,2),columnar_scores(:,2)) / (norm(adh_scores(:,2)) * norm(columnar_scores(:,2))));
variance_matrix(3,4) = acosd(dot(fea_scores(:,2),columnar_scores(:,2)) / (norm(fea_scores(:,2)) * norm(columnar_scores(:,2))));
%}

%hausdorff distance
%{
variance_matrix(1,1) = HausdorffDist(normal_scores(:,2),ones(size(normal_scores(:,2)))*mean(normal_scores(:,2)));
variance_matrix(2,2) = HausdorffDist(adh_scores(:,2),ones(size(adh_scores(:,2)))*mean(adh_scores(:,2)));
variance_matrix(3,3) = HausdorffDist(fea_scores(:,2),ones(size(fea_scores(:,2)))*mean(fea_scores(:,2)));
variance_matrix(4,4) = HausdorffDist(columnar_scores(:,2),ones(size(columnar_scores(:,2)))*mean(columnar_scores(:,2)));

variance_matrix(1,2) = HausdorffDist(normal_scores(:,2), adh_scores(:,2));
variance_matrix(1,3) = HausdorffDist(normal_scores(:,2), fea_scores(:,2));
variance_matrix(1,4) = HausdorffDist(normal_scores(:,2), columnar_scores(:,2));
variance_matrix(2,3) = HausdorffDist(adh_scores(:,2), fea_scores(:,2));
variance_matrix(2,4) = HausdorffDist(adh_scores(:,2), columnar_scores(:,2));
variance_matrix(3,3) = HausdorffDist(fea_scores(:,2), columnar_scores(:,2));
%}
t = num2cell(variance_matrix);
t = cellfun(@num2str, t, 'UniformOutput', false); 
X = [1 2 3 4 2 3 4 3 4 4];
Y = [1 1 1 1 2 2 2 3 3 4];
t = [t(1,1) t(1,2) t(1,3) t(1,4) t(2,2) t(2,3) t(2,4) t(3,3) t(3,4) t(4,4)]; 
figure; imagesc(variance_matrix); text(X, Y, t, 'HorizontalAlignment', 'Center'); colorbar;
title(sprintf('Taper Symmetry Class Variances\n Normal               Columnar               FEA                 ADH'));
ylabel('ADH             FEA             Columnar               Normal');

%separation
%euclidian distance
variance_matrix = zeros(4);

variance_matrix(1, 1) = norm(normal_scores(:, 3) - (ones(size(normal_scores(:,3)))*mean(normal_scores(:,3))));
variance_matrix(2, 2) = norm(columnar_scores(:, 3) - (ones(size(columnar_scores(:,3)))*mean(columnar_scores(:,3))));
variance_matrix(3, 3) = norm(fea_scores(:, 3) - (ones(size(fea_scores(:,3)))*mean(fea_scores(:,3))));
variance_matrix(4, 4) = norm(adh_scores(:, 3) - (ones(size(adh_scores(:,3)))*mean(adh_scores(:,3))));

variance_matrix(1, 2) = norm(normal_scores(:, 3) - columnar_scores(:, 3));
variance_matrix(1, 3) = norm(normal_scores(:, 3) - fea_scores(:, 3));
variance_matrix(1, 4) = norm(normal_scores(:, 3) - adh_scores(:, 3));
variance_matrix(2, 3) = norm(columnar_scores(:, 3) - fea_scores(:, 3));
variance_matrix(2, 4) = norm(columnar_scores(:, 3) - adh_scores(:, 3));
variance_matrix(3, 4) = norm(fea_scores(:, 3) - adh_scores(:, 3));


%vector angle
%{
variance_matrix(1, 1) = acosd(dot(normal_scores(:,3),ones(size(normal_scores(:,3)))*mean(normal_scores(:,3)))/(norm(normal_scores(:,3))*norm(ones(size(normal_scores(:,3)))*mean(normal_scores(:,3)))));
variance_matrix(2, 2) = acosd(dot(adh_scores(:,3),ones(size(adh_scores(:,3)))*mean(adh_scores(:,3)))/(norm(adh_scores(:,3))*norm(ones(size(adh_scores(:,3)))*mean(adh_scores(:,3)))));
variance_matrix(3, 3) = acosd(dot(fea_scores(:,3),ones(size(fea_scores(:,3)))*mean(fea_scores(:,3)))/(norm(fea_scores(:,3))*norm(ones(size(fea_scores(:,3)))*mean(fea_scores(:,3)))));
variance_matrix(1, 1) = acosd(dot(columnar_scores(:,3),ones(size(columnar_scores(:,3)))*mean(columnar_scores(:,3)))/(norm(columnar_scores(:,3))*norm(ones(size(columnar_scores(:,3)))*mean(columnar_scores(:,3)))));

variance_matrix(1,2) = acosd(dot(normal_scores(:,3),adh_scores(:,3)) / (norm(normal_scores(:,3)) * norm(adh_scores(:,3))));
variance_matrix(1,3) = acosd(dot(normal_scores(:,3),fea_scores(:,3)) / (norm(normal_scores(:,3)) * norm(fea_scores(:,3))));
variance_matrix(1,4) = acosd(dot(normal_scores(:,3),columnar_scores(:,3)) / (norm(normal_scores(:,3)) * norm(columnar_scores(:,3))));
variance_matrix(2,3) = acosd(dot(adh_scores(:,3),fea_scores(:,3)) / (norm(adh_scores(:,3)) * norm(fea_scores(:,3))));
variance_matrix(2,4) = acosd(dot(adh_scores(:,3),columnar_scores(:,3)) / (norm(adh_scores(:,3)) * norm(columnar_scores(:,3))));
variance_matrix(3,4) = acosd(dot(fea_scores(:,3),columnar_scores(:,3)) / (norm(fea_scores(:,3)) * norm(columnar_scores(:,3))));
%}

%hausdorff distance
%{
variance_matrix(1,1) = HausdorffDist(normal_scores(:,3),ones(size(normal_scores(:,3)))*mean(normal_scores(:,3)));
variance_matrix(2,2) = HausdorffDist(adh_scores(:,3),ones(size(adh_scores(:,3)))*mean(adh_scores(:,3)));
variance_matrix(3,3) = HausdorffDist(fea_scores(:,3),ones(size(fea_scores(:,3)))*mean(fea_scores(:,3)));
variance_matrix(4,4) = HausdorffDist(columnar_scores(:,3),ones(size(columnar_scores(:,3)))*mean(columnar_scores(:,3)));

variance_matrix(1,2) = HausdorffDist(normal_scores(:,3), adh_scores(:,3));
variance_matrix(1,3) = HausdorffDist(normal_scores(:,3), fea_scores(:,3));
variance_matrix(1,4) = HausdorffDist(normal_scores(:,3), columnar_scores(:,3));
variance_matrix(2,3) = HausdorffDist(adh_scores(:,3), fea_scores(:,3));
variance_matrix(2,4) = HausdorffDist(adh_scores(:,3), columnar_scores(:,3));
variance_matrix(3,3) = HausdorffDist(fea_scores(:,3), columnar_scores(:,3));
%}
t = num2cell(variance_matrix);
t = cellfun(@num2str, t, 'UniformOutput', false); 
X = [1 2 3 4 2 3 4 3 4 4];
Y = [1 1 1 1 2 2 2 3 3 4];
t = [t(1,1) t(1,2) t(1,3) t(1,4) t(2,2) t(2,3) t(2,4) t(3,3) t(3,4) t(4,4)]; 
figure; imagesc(variance_matrix); text(X, Y, t, 'HorizontalAlignment', 'Center'); colorbar;
title(sprintf('Separation Class Variances\n Normal               Columnar               FEA                 ADH'));
ylabel('ADH             FEA             Columnar               Normal');

%rows of X are the 24 observations (images), columns are top num_nuclei
%average ribbon, taper, and separation scores 
%{
X = [normal_scores ; adh_scores ; fea_scores ; columnar_scores];
Z = zscore(X);
[coefs, scores] = pca(Z);
vbls = {'Ribbon','Taper','Separation'}; % Labels for the variables
figure; biplot(coefs(:,1:2),'Scores',scores(1:6,1:2),'VarLabels',vbls,'MarkerEdgeColor','r','MarkerSize',12,'DisplayName','Normal');
hold on; biplot(coefs(:,1:2),'Scores',scores(7:12,1:2),'MarkerEdgeColor','k','MarkerSize',12);
hold on; biplot(coefs(:,1:2),'Scores',scores(13:18,1:2),'MarkerEdgeColor','g','MarkerSize',12);
hold on; biplot(coefs(:,1:2),'Scores',scores(19:24,1:2),'MarkerEdgeColor','m','MarkerSize',12);
%hack for the legend
h = zeros(4, 1);
h(1) = plot(NaN,NaN,'or');
h(2) = plot(NaN,NaN,'ok');
h(3) = plot(NaN,NaN,'og');
h(4) = plot(NaN,NaN,'om');
legend(h, 'Normal','ADH','FEA','Columnar');
title('Average scores for top 100 scored nuclei');
%}
%save('euc_variance_k12_100nuc_3lumen');

%lumen variances
%{
variance_matrix = zeros(4);

variance_matrix(1, 1) = norm(normal_scores(:, 4) - (ones(size(normal_scores(:,4)))*mean(normal_scores(:,4))));
variance_matrix(2, 2) = norm(adh_scores(:, 4) - (ones(size(adh_scores(:,4)))*mean(adh_scores(:,4))));
variance_matrix(3, 3) = norm(fea_scores(:, 4) - (ones(size(fea_scores(:,4)))*mean(fea_scores(:,4))));
variance_matrix(4, 4) = norm(columnar_scores(:, 4) - (ones(size(columnar_scores(:,4)))*mean(columnar_scores(:,4))));

variance_matrix(1, 2) = norm(normal_scores(:, 4) - adh_scores(:, 4));
variance_matrix(1, 3) = norm(normal_scores(:, 4) - fea_scores(:, 4));
variance_matrix(1, 4) = norm(normal_scores(:, 4) - columnar_scores(:, 4));
variance_matrix(2, 3) = norm(adh_scores(:, 4) - fea_scores(:, 4));
variance_matrix(2, 4) = norm(adh_scores(:, 4) - columnar_scores(:, 4));
variance_matrix(3, 4) = norm(fea_scores(:, 4) - columnar_scores(:, 4));

t = num2cell(variance_matrix);
t = cellfun(@num2str, t, 'UniformOutput', false); 
X = [1 2 3 4 2 3 4 3 4 4];
Y = [1 1 1 1 2 2 2 3 3 4];
t = [t(1,1) t(1,2) t(1,3) t(1,4) t(2,2) t(2,3) t(2,4) t(3,3) t(3,4) t(4,4)]; 
figure; imagesc(variance_matrix); text(X, Y, t, 'HorizontalAlignment', 'Center'); colorbar;
title(sprintf('Ribbon symmetry within/between class var LUMEN\n Normal               ADH               FEA                 Columnar'));
ylabel('Columnar             FEA             ADH               Normal');
%}
