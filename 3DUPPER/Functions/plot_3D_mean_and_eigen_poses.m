function plot_3D_mean_and_eigen_poses(Data_3D_align, Mean_pPCA, eignValues, eignVectors)

% define size variables
[Np,Framedim,Ns] = size(Data_3D_align);

%reshape Mean_pPCA into correct form (body-markers x x,y,z-coordinates)
mean_pose_3D_ppca = reshape(Mean_pPCA,[Np,Framedim]);


% FIGURE 1: plotting all selected frames and and mean pose
    % Define colors for each body marker
    colors = [...
        0.5, 0, 0.5;        % Nose - purple
        0, 0, 0.5;          % Left ear - dark blue
        0.5, 0.8, 1;        % Right ear - light blue
        0, 0.5, 0;          % Neck base - green
        0.5, 1, 0.5;        % Body anterior - light green
        1, 1, 0;            % Body posterior - yellow
        1, 0.5, 0];         % Tail base - orange
    
    % Create figure
    figure; hold on;
    
    % Loop over each body marker
    for markerIdx = 1:7
        % Extract all 4500 x,y,z positions for this marker
        marker_data = squeeze(Data_3D_align(markerIdx, :, :))';  % [4500 x 3]
    
        % Plot with scatter3, small size, 50% transparency
        scatter3(marker_data(:,1), marker_data(:,2), marker_data(:,3), ...
            15, colors(markerIdx,:), 'filled', 'MarkerFaceAlpha', 0.1);
    end
    
    % Plot mean coordinates as larger black dots
    scatter3(mean_pose_3D_ppca(:,1), mean_pose_3D_ppca(:,2), mean_pose_3D_ppca(:,3), ...
        60, 'k', 'filled');
    
    % Plot settings
    xlabel('X'); ylabel('Y'); zlabel('Z');
    grid on;
    view(3);
    xlim([-30 0]); ylim([0 30]); zlim([-10 10])
        
    title('3D Mouse Body Markers with Mean Pose');
            
    legend_labels = {'Nose', 'Left Ear', 'Right Ear', 'Neck Base', 'Body Anterior', 'Body Posterior', 'Tail Base', 'Mean Pose'};
    % Dummy handles for legend
    h = zeros(8,1);
    for i = 1:7
        h(i) = scatter3(NaN,NaN,NaN, 30, colors(i,:), 'filled');  % invisible dummy for legend
    end
    h(8) = scatter3(NaN,NaN,NaN, 60, 'k', 'filled');
    legend(h, legend_labels, 'Location', 'best');

    hold off



% FIGURE 2: visualizing eigenposes in animation
% Parameters
T = 100;               % Number of animation frames per eigenpose
Nbp = 7;               % Number of body parts (markers)
scale_factor = 3;      % Scaling factor for visualizing eigenvectors

% Reshape mean pose to a column vector [21 x 1]
mean_pose_vec = reshape(mean_pose_3D_ppca', [], 1);  % [21 x 1], column vector

% Create figure
figure;
h = subplot(1,1,1); hold on;

for n = 1:3  % Visualize the first 3 eigenposes
    for m = 1:T
        % Compute pose variation along eigenvector n with sine modulation
        pose_vec = mean_pose_vec + ...
            scale_factor * sqrt(eignValues(n)) * eignVectors(:,n) * sin(2*pi*m/(0.5*T));
        
        % Reshape back to [7 x 3] for 3D coordinates
        pose_3D = reshape(pose_vec, [3, Nbp])';  % Transpose to [7 x 3]
        
        % Plot 3D pose
        scatter3(pose_3D(:,1), pose_3D(:,2), pose_3D(:,3), ...
            60, 'r', 'filled');
        
        % Set plot limits and labels
         xlim([-20 -5]); ylim([5 20]); zlim([-5 5]);
        xlabel('X'); ylabel('Y'); zlabel('Z');
        title(['Eigenpose ' num2str(n)]);
        view(3);
        grid on;
        
        drawnow;
        pause(0.05);
        cla;  % Clear axes for next frame
    end
end
hold off


end