function[] = plot3d_video(X)
[Np,~,Ntt] = size(X);
figure; 
for t = 1:Ntt
    hold on;
    for m = 1:Np
        plot3(squeeze(X(m,1,t)),squeeze(X(m,2,t)),squeeze(X(m,3,t)),'.','MarkerSize',20);
    end
    xlim([-25 25]); ylim([-25 25]); zlim([-5 25]);
    view([-20 30]); title(num2str(t));
    %view([0 90]); 
    drawnow; 
    pause(0.04); cla;
end