function[] = Estimate_SSM_3D_(which_part, triangulated_data_path, SSM_model_name_path)
%INPUT: which_part is 'body' or 'tail' depending whether you re going to
%estimate Statistical Shape Model for body or tail. They will be saved separately. 

%initial parameters
Npose = 4500; %number of poses to estimate SSM, 4000-5000 is appropriate

% extact all filenames of triangulated data
all_triangulated_files = dir(fullfile(triangulated_data_path, '*mouse*'));
filename = cell(1,length(all_triangulated_files));

% loop through all filenames and load them into cell array
for i = 1:length(all_triangulated_files)
    
    filename{i} = all_triangulated_files(i).name;

end

%MIKE: Input body part
if strcmp(which_part,'body')
   ind = [1:7];
elseif strcmp(which_part,'tail')
   ind = [7:10]; 
else
   disp('Valid INPUTs: "body" or "tail"');
   return
end

%load 3D coordinates from all files
Nfile = numel(filename);
Maindata = [];
for n = 1:Nfile
    
    n_data = load(fullfile(triangulated_data_path,filename{n}), 'X');
    n_data = n_data.X;

    Maindata = cat(3,Maindata, n_data);
 
end

selected_3D_data = Maindata;
clear Maindata;


%remove some body parts (body or tail)
selected_3D_data = selected_3D_data(ind,:,:);
Nbp = numel(ind);

%To generate the model identify poses with no outliers and select a random
%selction of Npose from them;

ind_rand = randperm(size(selected_3D_data,3), Npose);
selected_3D_data = selected_3D_data(:,:,ind_rand);

%%%
% Estimation_Model function eceives the 3D data and return the Data which all missing data filled up and  mean
% estimated by Ransac, Mean of pPCA, Covariance of pPCA and Eigenvalues
% and eigenpose
[Data_3D_KNN Mean_Ransac_3D Mean_pPCA Cov_pPCA eignValues eignVectors]=Estimation_Model(selected_3D_data,0.8);


% save model
save(SSM_model_name_path, "Data_3D_KNN", "Mean_Ransac_3D", "Mean_pPCA", "Cov_pPCA", "eignValues", "eignVectors"); 
