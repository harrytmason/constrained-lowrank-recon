function [T,rT,nT, fT,rTvec] = CLRtsub(D,E,X,T,lamt,Tp,lams,tol,Nin,Nt,r)
%% Calcuates the solution to the temporal subproblem of k-t PERRI
%
% Inputs:  D:   nufft_adj(k-space data)
%          E:   nufft structure + more 
%          X:   current guess of X
%          T:   current guess of T
%          lamt: the weighting on the temporal constraint
%          lams: the weighting on the smoothness constraint
%          Tp:  the prior of T
%          tol: tolerance of problem
%          Nin: number of iterations of subproblem  
%
% Outputs: T:     optimised T
%          rT:    residual output of cost function
%          nT:    number of iterations of subproblem 
%          fT:    state of subproblem upon completion
%          rTVec: vector of residual outputs
%
%%
%% Code
             TFun  = @(x)vec(X'*E.mtimes2(X*reshape(x,[r, Nt]))) + lamt*x + lams * vec([zeros(r,1) diff(reshape(x,[r, Nt]),1,2) ]+[  fliplr(diff(fliplr(reshape(x,[r, Nt])),1,2)) zeros(r,1)]);
[T,fT,rT,nT,rTvec] = minres(TFun, vec(X'*D + lamt*Tp'), tol, Nin, [], [], vec(T'));
             T     = reshape(T, [r, Nt])';
end