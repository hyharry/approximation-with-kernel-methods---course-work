function visualize(f, Xvis, Yvis, Xcenter, gamma)
% f: approx function
% Xvis: visualization points
% Xcenter: center for f evaluation
% gamma: kernel parameter
% plot_prec: plot parameter

% $Author: yihu $	$Date: 2016/01/21 15:00:30 $	$Revision: 0.1 $

% x = Xvis(:,1);
% y = Xvis(:,2);
% z = f(Xvis);
% d = plot_prec;
% 
% % Generate griddata using vector x,y
% [Xmesh,Ymesh]=meshgrid(min(x):d:max(x),min(x):d:max(x));
% 
% row_ind = Xvis(:,1)>0 & Xvis(:,2)>0;
% Xvis(row_ind,:) = [];

% Zmesh = griddata(x,y,z,Xmesh,Ymesh);

z = f([Xvis(:), Yvis(:)]);
Zvis = reshape(z,size(Xvis));

% Plot
figure
surf(Xvis,Yvis,Zvis);
% surf(Xmesh,Ymesh,Zmesh);
% surfc(Xmesh,Ymesh,Zmesh);
% surf(Xmesh,Ymesh,Zmesh);
% contourf(Xmesh,Ymesh,Zmesh);
hold on
% scatter3(Xcenter(:,1), Xcenter(:,2),f(Xcenter),'r','filled');
scatter(Xcenter(:,1), Xcenter(:,2),'r','filled');

if gamma == 0
    title('exact');
else
    title(['\gamma =', num2str(gamma)]);
end

end