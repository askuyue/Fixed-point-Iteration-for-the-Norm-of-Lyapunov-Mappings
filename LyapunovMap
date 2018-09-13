
clear
clc
format long
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
        max_{X \neq 0} || AX + XA' ||_{F}^{2} / || X ||_{F}^2


%}

% 矩阵的阶数 N
N                           =               6;
A                           =               [1,2,0,2,3,0;
                                                0,3,0,2,0,1;
                                                2,2,3,0,0,3;
                                                1,2,0,3,1,0;
                                                2,0,0,0,1,3;
                                                0,1,0,0,2,1];
                                        
% A           =                   hilb(N);

%  A   =   randn(N,N);


%A           =           A ./ norm(A,'fro');

% 单位矩阵
IN                          =               eye(N);

% Lyapunov Matrix : LA 
LA                          =               kron(A,IN) + kron(IN,A);

% 求解矩阵LA的最大奇异值及其对应的奇异向量
[V,D]               = eig(LA'*LA);
    

% 找到LA'*LA的最大特征值及其特征向量
diagD  = diag(D);
[LambdaMax,Idmax]           =           max(diagD);
VectorV                              =           V(:,Idmax);
MatrixV                              =           reshape(VectorV',[N,N]);

% 测试
fprintf("The Maximum Eigenvalues of LA'LA : %f\n",LambdaMax);
LMAx = norm(A*MatrixV + MatrixV*A','fro')^2 / norm(MatrixV,'fro')^2;
fprintf("The Maximum Eigenvalues of LA'LA : %f\n",LMAx);


IterMax                 =       10;

Alambda               =        A ./ sqrt(LambdaMax);
% 每一步产生的迭代矩阵的F范数
FnormIters   =   zeros(IterMax,1);

% 每一步产生的迭代矩阵的目标函数值
ObjIters        =   zeros(IterMax,1);

%  迭代收敛速率的考虑
% % || X^{k+1} - Xs ||_{F} / || X^{k} - Xs ||_{F}
% RateIters        =   zeros(IterMax,1);

% 目标函数收敛速率的考虑
% (F* - FObjIters) / F*
ObjRateIters        =   zeros(IterMax,1);

% 把迭代序列单位化之后与解之间的关系
XmatrixK        =   zeros(IterMax,1);

% 把迭代序列单位化之后与解之间的关系
dXmatrixK        =   zeros(IterMax,1);
% FM = FixMap(X,A)
% % % 初始矩阵
%     X0              =   rand(N,N) ./ N ;
% % 对称化
%     X0              =   (X0 + X0') ./ 2;
    
    X0  =   hilb(N);
    X0  =   X0 ./ norm(X0,'fro');
    
    Xk              =   X0;
    Xkp1          =   Xk;
    
    k = 1;
while k < (IterMax+1)
     
    Xmid         =  Xkp1;
    % 把迭代序列单位化之后与解之间的关系
    XmatrixK(k)        =   norm(Xmid ./ norm(Xmid,'fro') - MatrixV, 'fro') ;
    dXmatrixK(k)      =   norm(Xmid ./ norm(Xmid,'fro') + MatrixV, 'fro') ;
    % 每一步产生的迭代矩阵的F范数
    FnormIters(k)    =   norm(Xmid,'fro');
    % 每一步产生的迭代矩阵的目标函数值
    ObjIters(k)        =   norm(A*Xmid + Xmid*A','fro')^2 / norm(Xmid,'fro')^2;
    % 目标函数收敛速率的考虑
    % (F* - FObjIters) / F*
    ObjRateIters(k)        =   (LambdaMax - ObjIters(k)) / LambdaMax;
    % 计算下一步的迭代矩阵
    Xkp1         =  FixMap(Xmid,Alambda);
% %     % 每一步产生的迭代矩阵的F范数
% %     FnormIters(k)    =   norm(Xkp1,'fro');
% %     % 每一步产生的迭代矩阵的目标函数值
% %     ObjIters(k)        =   norm(A*Xkp1 + Xkp1*A','fro')^2 / norm(Xkp1,'fro')^2;
% %     % 目标函数收敛速率的考虑
% %     % (F* - FObjIters) / F*
% %     ObjRateIters(k)        =   (LambdaMax - ObjIters(k)) / LambdaMax;

%     %  迭代收敛速率的考虑
%     % || X^{k+1} - Xs ||_{F} / || X^{k} - Xs ||_{F}
%     alpha  = 1;
%     RateIters(k)        =   norm(Xkp1 ./ norm(Xkp1,'fro') - MatrixV,'fro') / norm(Xmid ./ norm(Xmid,'fro') - MatrixV,'fro')^(alpha);
%    
        k = k+1;
end

Xkp1            =           Xkp1 ./ norm(Xkp1,'fro');

% 每一步产生的迭代矩阵的F范数
figure(1)
plot(1:IterMax,FnormIters,'b:d','LineWidth',1.5)
xlabel('Iterations')
ylabel('|| X^{k} ||_{F}')
axis([1 IterMax min(FnormIters)-0.01 max(FnormIters)]);

% 每一步产生的迭代矩阵的目标函数值
figure(2)
plot(1:IterMax,ObjIters,'g--h','LineWidth',1.5)
hold on
plot(1:IterMax,LambdaMax.*ones(IterMax,1),'r--h','LineWidth',1.5);
xlabel('Iterations')
ylabel('g(X^{k})')
legend('g(X^{k})','\lambda_{max}')
axis([1 IterMax min(ObjIters) max(ObjIters)+2]);

% ff = zeros(IterMax,1);
% for j = 1:IterMax
%     ff(j) = 1 / j^5;
% end
 % 目标函数收敛速率的考虑
 % (F* - FObjIters) / F*
%  figure(3)
%  plot(1:IterMax,ObjRateIters,'m-s','LineWidth',1.5)
% %  hold on
% %  plot(1:IterMax,ff,'c-.+')
%  xlabel('Iterations')
%  ylabel('(\lambda_{max} - g(X^{k}) )/ \lambda_{max}')
%  axis([1 IterMax min(ObjRateIters)-0.01 max(ObjRateIters)]);
% %  
% %  迭代收敛速率的考虑
% % || X^{k+1} - Xs ||_{F} / || X^{k} - Xs ||_{F}
% figure(4)
% plot(1:IterMax,RateIters,'k:p','LineWidth',2)
% xlabel('Iterations')
% ylabel('|| X^{k+1} / || X^{k+1} ||_{F}- X^{*} || / || X^{k} / || X^{k} ||_{F} - X^{*} ||')


% % 把迭代序列单位化之后与解之间的关系
figure(5)
plot(1:IterMax,XmatrixK,'k-.o','LineWidth',1.5)
xlabel('Iterations')
ylabel('|| X^{k}  - X^{*} ||_{F} / || X^{*} ||_{F} ')
axis([1 IterMax min(XmatrixK) max(XmatrixK)]);

% % 把迭代序列单位化之后与解之间的关系
figure(6)
plot(1:IterMax,dXmatrixK,'m:o','LineWidth',2)
xlabel('Iterations')
ylabel('|| X^{k} + X^{*} ||_{F} / || X^{*} ||_{F} ')
axis([1 IterMax min(dXmatrixK) 2]); %max(dXmatrixK)



