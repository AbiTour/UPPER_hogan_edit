function[Xfit,b,A,T,C,missing] = fit_data_x(X,lambda,mean_pose,P,var_res,min_num)
alpha_reg = 0.1;% 0.001;
Nshape = numel(lambda);
Nbp = size(X,1);
stop_search = false;
options = optimoptions('fminunc','Display','none');
%init shape parameters
b0 = zeros(Nshape,1);
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