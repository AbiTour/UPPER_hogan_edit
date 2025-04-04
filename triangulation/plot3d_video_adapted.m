function[] = plot3d_video_adapted(Data_3D_KNN, make_lines)

[Np,~,Ntt] = size(Data_3D_KNN);
figure; 
for t = 1:1:100
    hold on;
    for m = 1:Np
        plot3(squeeze(Data_3D_KNN(m,1,t)),squeeze(Data_3D_KNN(m,2,t)),squeeze(Data_3D_KNN(m,3,t)),'.','MarkerSize',20);
    end
    %
    if make_lines
        line([Data_3D_KNN(1,1,t) Data_3D_KNN(3,1,t)],[Data_3D_KNN(1,2,t) Data_3D_KNN(3,2,t)],[Data_3D_KNN(1,3,t) Data_3D_KNN(3,3,t)],'Color','k','LineWidth',2);
        line([Data_3D_KNN(2,1,t) Data_3D_KNN(3,1,t)],[Data_3D_KNN(2,2,t) Data_3D_KNN(3,2,t)],[Data_3D_KNN(2,3,t) Data_3D_KNN(3,3,t)],'Color','k','LineWidth',2);
        %
        line([Data_3D_KNN(3,1,t) Data_3D_KNN(7,1,t)],[Data_3D_KNN(3,2,t) Data_3D_KNN(3,2,t)],[Data_3D_KNN(3,3,t) Data_3D_KNN(7,3,t)],'Color','k','LineWidth',2);
        line([Data_3D_KNN(3,1,t) Data_3D_KNN(8,1,t)],[Data_3D_KNN(3,2,t) Data_3D_KNN(8,2,t)],[Data_3D_KNN(3,3,t) Data_3D_KNN(8,3,t)],'Color','k','LineWidth',2);
        line([Data_3D_KNN(7,1,t) Data_3D_KNN(9,1,t)],[Data_3D_KNN(7,2,t) Data_3D_KNN(9,2,t)],[Data_3D_KNN(7,3,t) Data_3D_KNN(9,3,t)],'Color','k','LineWidth',2);
        line([Data_3D_KNN(8,1,t) Data_3D_KNN(10,1,t)],[Data_3D_KNN(8,2,t) Data_3D_KNN(10,2,t)],[Data_3D_KNN(8,3,t) Data_3D_KNN(10,3,t)],'Color','k','LineWidth',2);
        %
        line([Data_3D_KNN(7,1,t) Data_3D_KNN(8,1,t)],[Data_3D_KNN(7,2,t) Data_3D_KNN(8,2,t)],[Data_3D_KNN(7,3,t) Data_3D_KNN(8,3,t)],'Color','k','LineWidth',2);
        line([Data_3D_KNN(9,1,t) Data_3D_KNN(10,1,t)],[Data_3D_KNN(9,2,t) Data_3D_KNN(10,2,t)],[Data_3D_KNN(9,3,t) Data_3D_KNN(10,3,t)],'Color','k','LineWidth',2);
        line([Data_3D_KNN(9,1,t) Data_3D_KNN(11,1,t)],[Data_3D_KNN(9,2,t) Data_3D_KNN(11,2,t)],[Data_3D_KNN(9,3,t) Data_3D_KNN(11,3,t)],'Color','k','LineWidth',2);
        line([Data_3D_KNN(10,1,t) Data_3D_KNN(12,1,t)],[Data_3D_KNN(10,2,t) Data_3D_KNN(12,2,t)],[Data_3D_KNN(10,3,t) Data_3D_KNN(12,3,t)],'Color','k','LineWidth',2);
        line([Data_3D_KNN(11,1,t) Data_3D_KNN(12,1,t)],[Data_3D_KNN(11,2,t) Data_3D_KNN(12,2,t)],[Data_3D_KNN(11,3,t) Data_3D_KNN(12,3,t)],'Color','k','LineWidth',2);
        %
        line([Data_3D_KNN(3,1,t) Data_3D_KNN(4,1,t)],[Data_3D_KNN(3,2,t) Data_3D_KNN(4,2,t)],[Data_3D_KNN(3,3,t) Data_3D_KNN(4,3,t)],'Color','k','LineWidth',2);
        line([Data_3D_KNN(4,1,t) Data_3D_KNN(5,1,t)],[Data_3D_KNN(4,2,t) Data_3D_KNN(5,2,t)],[Data_3D_KNN(4,3,t) Data_3D_KNN(5,3,t)],'Color','k','LineWidth',2);
        line([Data_3D_KNN(5,1,t) Data_3D_KNN(6,1,t)],[Data_3D_KNN(5,2,t) Data_3D_KNN(6,2,t)],[Data_3D_KNN(5,3,t) Data_3D_KNN(6,3,t)],'Color','k','LineWidth',2);
    end
    %
    if make_lines
        xlim(0.6*[-15 15]); ylim(0.6*[-15 15]); zlim(0.6*[-15 15]);
    else
        xlim([-25 25]); ylim([-25 25]); zlim([-5 20]);
    end
    view([-70 20]); title(num2str(t));
    % view([45 0]); 
    % xlim([10 20])
    % ylim([10 20])
    % zlim([0 5])
     xlim auto
     ylim auto
     zlim auto
    drawnow; 
    pause(0.4); 
    if t ~= Ntt
        cla;
    end
end