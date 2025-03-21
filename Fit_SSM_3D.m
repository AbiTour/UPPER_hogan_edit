function[Xfit,b,t] = Fit_SSM_3D(triangulated_data_path, which_part, save_path)

%which_part is 'body' or 'tail' depending whether you re going to
%use fit SSM for body or tail. They will be saved separately. 

%init parameters
Ndim = 3; 
THeig = 0.9;


all_triangulated_files = dir(fullfile(triangulated_data_path, '*mouse*'));
filename = cell(1,length(all_triangulated_files));

% loop through all filenames and load them into cell array
for i = 1:length(all_triangulated_files)
    
    filename{i} = all_triangulated_files(i).name;
end

% Load SSM
if strcmp(which_part,'body')
   ind = 1:7;
   load('mouse1_10_3Dssm.mat','Mean_pPCA','eignValues','eignVectors');
   min_num = 6;
elseif strcmp(which_part,'tail')
   ind = 8:10; 
   load('SSM_tail.mat','Mean_pPCA','eignValues','eignVectors');
   min_num = 3;
else
   disp('Valid INPUTs: "body" or "tail"');
   return
end

%determine the number of eigenvalues
Neig = min(find(cumsum(eignValues)/sum(eignValues)>THeig));

%load coordinates
%load 3D coordinates from all files
Nfile = numel(filename);
Maindata = [];
for n = 1:Nfile
    
    n_data = load(fullfile(triangulated_data_path,filename{n}), 'X');
    n_data = n_data.X;

    Maindata = cat(3,Maindata, n_data);
 
end

Nframe = size(Maindata, 3);

%select only body or tail keypoints
X_all_data = Maindata(ind,:,:);
Nbp = numel(ind);

%rearrange model 
lambda = eignValues(1:Neig);
var_res = mean(eignValues(Neig+1:end));
mean_pose = reshape(Mean_pPCA, [Nbp, Ndim]); %mean_pose = Mean_pPCA;
eigen2=reshape(eignVectors(:,:),Ndim,Nbp,Nbp*Ndim);
for i=1:Nbp*Ndim
    P(:,:,i)=(eigen2(:,:,i)');
end
P = P(:,:,1:Neig);

%fit the 3D data
Xfit_all_data = zeros(Nbp,Ndim,Nframe); 
b_all_data = zeros(Neig,Nframe);
C_all_data = zeros(1,Nframe);
A_all_data = zeros(1,Nframe);
T_all_data = zeros(3,Nframe);
missing_all_data = true(1,Nframe);

parfor n = 1:Nframe
    [Xfit_all_data(:,:,n),b_all_data(:,n),A_all_data(n),T_all_data(:,n),C_all_data(n),missing_all_data(n)] = fit_data(X_all_data(:,:,n),lambda,mean_pose,P,var_res,min_num);
    fprintf('frame %s of %s\n',num2str(n),num2str(Nframe));
end

% Save data separately for each file
frame_start = 1;

for n = 1:Nfile
    % Load the original data to get the frame count for each file
    n_data = load(fullfile(triangulated_data_path, filename{n}), 'X');
    n_data = n_data.X;
    
    % Extract number of frames per file
    Nframe_file = size(n_data, 3);
    % Calculate final frame index of current file (with respect to
    % concatenated variables)
    frame_end = frame_start + Nframe_file - 1;

    % Extract corresponding portion of fitted data
    Xfit = Xfit_all_data(:, :, frame_start:frame_end);
    b = b_all_data(:, frame_start:frame_end);
    C = C_all_data(frame_start:frame_end);
    X = X_all_data(:, :, frame_start:frame_end);
    missing = missing_all_data(frame_start:frame_end);
    A = A_all_data(:, frame_start:frame_end);
    T = T_all_data(:, frame_start:frame_end);

    % Save the fitted data with the appropriate filename
    save(fullfile(save_path, [filename{n}(1:end-4) '_' which_part '_fit.mat']), ...
        'b', 'Xfit', 'C', 'X', 'missing', 'A', 'T');

    % Update frame counter
    frame_start = frame_end + 1;
end


function[Xfit,b,A,T,C,missing] = fit_data(X,lambda,mean_pose,P,var_res,min_num)
alpha_reg = 0.1;% 0.001;
Nshape = numel(lambda);
Nbp = size(X,1);
stop_search = false;
options = optimoptions('fminunc','Display','none');
%init shape parameters
b0 = zeros(Nshape,1);

% Check if there are any NaNs in X
    if sum(isnan(X(:)))>5
        missing = true;  % Mark this frame as missing
        Xfit = NaN * mean_pose;
        T = NaN(3, 1);  % Set 3D translation as NaN
        A = NaN;        % Set angle as NaN
        b = NaN * ones(Nshape, 1);
        C = 1000;       % High cost for missing data
        return;
    end
    
%fit the SSM
ind_num = find(~isnan(X(:,1)));
missing = false;
if numel(ind_num)>=min_num
    b = fminunc(@(b) fit_SSM(b, X(ind_num,:), mean_pose(ind_num,:), lambda, P(ind_num,:,:), alpha_reg, var_res), b0, options);
    [C, ~, R, T] = fit_SSM(b, X(ind_num,:), mean_pose(ind_num,:), lambda, P(ind_num,:,:), alpha_reg, var_res);
    Xfit = mean_pose;
    for n = 1:Nshape
        Xfit = Xfit + b(n)*P(:,:,n);
    end
    Xfit = Xfit*R + repmat(T,Nbp,1);
    A = rotm2eul(R);
    A = A(1);
else
    missing = true;
    Xfit = NaN*mean_pose;
    T = [NaN; NaN];
    A = NaN;
    b = NaN*ones(Nshape,1);
    C = 1000;
end

%show results
make_fig = false;
if make_fig & ~missing & ind_num<6
    fig1 = figure; hold on
    plot(X(:,1),X(:,2),'k.','MarkerSize',16); 
    plot(Xfit(:,1),Xfit(:,2),'bo','MarkerSize',8,'LineWidth',2); 
    title(['mean fit: C = ' num2str(C)]);
    mX = nanmean(X);
    xlim([mX(1)-150 mX(1)+150]);
    ylim([mX(2)-150 mX(2)+150]);
    ginput(); close all;
end

function[C, Xfit, R, T] = fit_SSM(b, X, mu, lambda, P, alpha_reg, var_res)
Nshape = numel(lambda);
%
Xfit0 = mu;
for n = 1:Nshape
    Xfit0 = Xfit0+b(n)*P(:,:,n);
end
%
[~, Xfit, tr] = procrustes(X,Xfit0,'Reflection',false, 'Scaling',false);
%
dist = sum((X(:)-Xfit(:)).^2)/var_res;
%
reg = alpha_reg*sum((b.^2)./lambda);
%
C = dist+reg;
%
R = tr.T;
T = tr.c(1,:);
%
