% Danny Drieﬂ
% Compares different ways to compute the squared distance of two matrices
% in terms of speed

clear;

nTest = [10^1,10^2,10^3,10^4];
timeBsxFun = zeros(size(nTest,2),1);
timeMatrix = zeros(size(nTest,2),1);
timePDist2 = zeros(size(nTest,2),1);
i = 1;
for n = nTest
    X = rand(2,n);
    XPrime = X;
    
    tic
    KBsxFun = (bsxfun(@plus,sum(X.^2,1)',sum(XPrime.^2,1)) - 2*X'*XPrime);
    t = toc;
    timeBsxFun(i,1) = t;
    
    tic
    KMatrix = sum(X.^2,1)'*ones(1,size(XPrime,2)) + ones(size(X,2),1)*sum(XPrime.^2,1) - 2*X'*XPrime;
    t = toc;
    timeMatrix(i,1) = t;
    
    tic 
    KPDist = pdist2(X',XPrime').^2;
    t = toc;
    timePDist2(i,1) = t;
    
    i = i+1; 
end;

loglog(nTest,timeBsxFun,nTest,timeMatrix,nTest,timePDist2)
legend('Bsxfun','Matrix','PDist2')
xlabel('n')
ylabel('time [s]')