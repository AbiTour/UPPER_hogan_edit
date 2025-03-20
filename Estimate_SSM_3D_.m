function[] = Estimate_SSM_3D_abi(which_part, all_DLC_files, SSM_model_name)
%INPUT: which_part is 'body' or 'tail' depending whether you re going to
%estimate Statistical Shape Model for body or tail. They will be saved separately. 

%initial parameters
TH = 0.9; %likelihood threshold for outliers
Npose = 4500; %number of poses to estimate SSM, 4000-5000 is appropriate
K = 3; %K nearest neighbour

% intialise cell array to contain all filenames
filename = cell(1,length(all_DLC_files));

% loop through all filenames and load them into cell array
for i = 1:length(all_DLC_files)
    
    filename{i} = all_DLC_files(i).name;

end

%MIKE: Input body part
if strcmp(which_part,'body')
   ind = [1:7 11 12];
elseif strcmp(which_part,'tail')
   ind = [7:10]; 
else
   disp('Valid INPUTs: "body" or "tail"');
   return
end

%load coordinates
[D_train,lik] = load_coordinates(filename);
[Nbp,Nframe] = size(lik);

%remove some body parts (body or tail)
D_train = D_train(ind,:,:);
lik = lik(ind,:); 
Nbp = numel(ind);

%To generate the model identify poses with no outliers and select a random
%selction of Npose from them;
ind_good = find(min(lik)>TH);
Ngood = numel(ind_good);
ind_rand = randperm(Ngood);
Npose = min(Npose,Ngood);
ind_rand = ind_rand(1:Npose);
D_train = D_train(1:7,:,ind_rand); % 1:7 

%%%
% Estimation_Model function eceives the 3D data and return the Data which all missing data filled up and  mean
% estimated by Ransac, Mean of pPCA, Covariance of pPCA and Eigenvalues
% and eigenpose
[Data_3D_KNN Mean_Ransac_3D Mean_pPCA Cov_pPCA eignValues eignVectors]=Estimation_Model( RawData3D_full,0.8);


%%

