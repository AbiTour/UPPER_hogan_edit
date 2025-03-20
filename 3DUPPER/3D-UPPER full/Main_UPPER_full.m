clear all; close all; clc; 
%% %%  Data reading  and preperation should be do by the user base on her/his data
load('D:\PhD 2nd Year\dlc_data_cohort2\triangulated_3D_data\threshold_06\mouse1_extinction_p1_3D_triangulated.mat')
RawData=X;
% [ND,Ns] = size(RawData);
% Framedim=3;
% Np=ND/Framedim;
% RawData3D=reshape(RawData,[Np,Framedim,Ns]); %reshape
RawData3D_full=RawData(1:7,:,3000:5000);
%%
%%%
% Estimation_Model function eceives the 3D data and return the Data which all missing data filled up and  mean
% estimated by Ransac, Mean of pPCA, Covariance of pPCA and Eigenvalues
% and eigenpose
[Data_3D_KNN, Mean_Ransac_3D, Mean_pPCA, Cov_pPCA, eignValues, eignVectors]=Estimation_Model( RawData3D_full,0.8);
%%
%the input is the 3D Rawdata (it is suggested use the Data_3D_KNN which has no missing values) and the output is 3D reconstructed data which
%backed to original place in the arena. 

[Reconstructed_Data_full]=Reconstruct_Data(RawData3D_full,Data_3D_KNN,0.99,Mean_Ransac_3D,Mean_pPCA,Cov_pPCA);

%%
save('Reconstructed_Data_full')


