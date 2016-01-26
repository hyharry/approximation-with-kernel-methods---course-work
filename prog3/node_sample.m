function [X, quality] = node_sample(n,nrand,a,b,c,d)
% Simple random sampling in a domain. Quality is not guaranteed

% X = [0 0;1 0;1 1;0 1];
x = linspace(0,1,nrand)';
X = [zeros(size(x)), x;
    x, zeros(size(x));
    ones(size(x)), x;
    x, ones(size(x))];
A = [1 0 1 0 1];

for i=1:n
    [~, ind] = max(A(:,1));
    Amax = A(ind,:);
    [Xadd, Asmall] = split(Amax);
    A(ind,:) = [];
    A = [A; Asmall];
    X = [X; Xadd];
end

quality = min(A(:,1))/max(A(:,1));

if nargin==6
    X(:,1) = a + (b-a)*X(:,1);
    X(:,2) = c + (d-c)*X(:,2);
end

X = unique(X,'rows');

tri = delaunay(X(:,1),X(:,2));
triplot(tri,X(:,1),X(:,2));

function [X, Asmall] = split(A)

a = A(2);
b = A(3);
c = A(4);
d = A(5);
x = a + (b-a)*rand;
y = c + (d-c)*rand;
X = [x y];
Asmall = [(x-a)*(y-c) a x c y;
          (b-x)*(d-y) x b c y;
          (x-a)*(d-y) a x y d;
          (b-x)*(d-y) x b y d];
