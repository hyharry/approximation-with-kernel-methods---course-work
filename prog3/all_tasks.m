function all_tasks

% Generate gamma
i = -1:3;
gamma_li = 2*10.^i;

% Finite quotients and direct derivative, 11 equal dist centers in [0,1]
Task2(10,gamma_li,'f');
Task2(10,gamma_li,'d');

% Finer centers but the same results
Task2(20,gamma_li,'f');
Task2(20,gamma_li,'d');

% alpha=0.5, Finer centers does not help for Task3
Task3(0.5,10,gamma_li,'f');
Task3(0.5,20,gamma_li,'f');
Task3(0.5,30,gamma_li,'f');
Task3(0.5,40,gamma_li,'f');
Task3(0.5,50,gamma_li,'f');

Task3(0.5,10,gamma_li,'d');
Task3(0.5,20,gamma_li,'d');
Task3(0.5,30,gamma_li,'d');
Task3(0.5,40,gamma_li,'d');

% alpha=0.5, error does not change, when gamma larger
Task3(0.5,20,[200,2000,20000,2e5,2e6],'f');

% Different alpha, alpha=1,2,4, direct method
Task3(1,30,gamma_li,'d');
Task3(2,30,gamma_li,'d');
Task3(4,20,gamma_li,'d'); % better approximation