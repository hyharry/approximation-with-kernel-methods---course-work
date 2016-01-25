function model = poisson_model(params)
%function model = poisson_model(params)
%
% small example of a model, i.e. a structure describing the data
% functions and geometry information of a general elliptic equation consisting 
% of diffusion, convection, reaction equation:
%
% - div ( a(x) grad u(x)) + div (b(x) u(x)) + c(x) u(x) = f(x)    on Omega
%                                                 u(x)) = g_D(x)  on Gamma_D
%                                       a(x) grad u(x)) = g_N(x)  on Gamma_N
% 
%  Here, we denote the functions as
%                   u: solution (if known, useful for validation purpose)
%                   f: source
%                   a: diffusivity
%                   b: velocity
%                   c: reaction
%                 g_D: Dirichlet boundary values
%                 g_N: Neumann boundary values
%
% Each function allows the evaluation in many points
% simultaneuously by
%
%        model.source(glob)
% or     model.source(glob,params)
%
% where glob is a n times 2 matrix of row-wise points. The result
% is a n times 1 vector of resulting values of f.
%
% Additionally, for some discretization methods, further functions
% are required, e.g. derivatives of the data functions:
%
%        model.diffusivity_gradient(glob,params)
%
% result is a matrix of size(glob), each row contains the gradient of a(x)
% in the point given in the corresponding row of glob.
%
% additionally, the model has a function, which determines, whether
% a point lies on a Dirichlet or Neumann boundary:
%        
%           model.boundary_type(glob)
%                0 no boundary (inner edge or point)
%               -1 indicates Dirichlet-boundary
%               -2 indicates Neumann-boundary
%
% Additionally, the normals in a boundary point can be requested by
%
%           model.normals(glob)
%           
% Here, glob are assumed to be boundary points lying ON THE
% INTERIOR of an edge, such that the outer unit normal is well-defined.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% The data functions given in this model are the simple poisson
% equation with Gamma_D = boundary(Omega), Gamma_N = {}
%
%    -div (grad u) = f
%            u = 0   on   Gamma_D
%
% with exact solution u(x) = x_1(1-x_1)x_2(1-x_2), 
% i.e. f(x) = 2 (x_1 + x_2 - x_1^2 - x_2^2)
%
% params is an optional parameter, perhaps useful later

% B. Haasdonk 23.11.2010

model = [];

% params is an optional parameter, perhaps useful later
model.source = @(glob,params) ...
            2 * (glob(:,1)+ glob(:,2)- glob(:,1).^2 - glob(:,2).^2);
model.reaction = @(glob,params) zeros(size(glob,1),1);
model.velocity = @(glob,params) zeros(size(glob,1),1);
model.diffusivity = @(glob,params) ones(size(glob,1),1);
model.diffusivity_gradient = @(glob,params) zeros(size(glob,1),2);

model.reaction = @(glob,params) zeros(size(glob,1),1);

%  Dirichlet everywhere, see function below.
model.boundary_type = @my_boundary_type;
model.dirichlet_values = @(glob,params) zeros(size(glob,1),1);
model.neumann_values = @(glob,params) zeros(size(glob,1),1);

model.normals = @my_normals;

% solution is known:
model.solution = @(glob,params) ...
             glob(:,1).*(1-glob(:,1)).*glob(:,2).*(1-glob(:,2));

% all edges of unit square are dirichlet, other inner
function res = my_boundary_type(glob,params)
res = zeros(size(glob,1),1);
i  = find(glob(:,1)<=1e-10);
i  = [i; find(glob(:,1)>=1-1e-10)];
i  = [i; find(glob(:,2)<=1e-10)];
i  = [i; find(glob(:,2)>=1-1e-10)];
res(i) = -1;

function res = my_normals(glob,params)
res = zeros(size(glob,1),2); % each row one normal
i  = find(glob(:,1)>1-1e-10);
res(i,1)= 1.0;
i  = find(glob(:,1)<+1e-10);
res(i,1)= -1.0;
i  = find(glob(:,2)>1-1e-10);
res(i,2)= 1.0;
i  = find(glob(:,2)<+1e-10);
res(i,2)= -1.0;

