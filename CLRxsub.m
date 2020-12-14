function [X,rX,nX, fX,rXvec] = CLRxsub(D,E,X,v,lamx,Xp, tol, Nin, Nx,r)
%% Calcuates the solution to the u subproblem of k-t PERRI
%
% Inputs:  D:    nufft_adj(k-space data)
%          E:    nufft structure + more 
%          X:    current guess of X
%          T:    current guess of T
%          lamx: the weighting on the spatial constraint
%          Xp:   the prior of X
%          tol:  tolerance of problem
%          Nin:  number of iterations of subproblem  
%          Nx:   number of pixels/volumes (dim1 of X)
%          r:    number of components (dim2 of X)
%
% Outputs: X:     optimised X
%          rX:    residual output
%          nX:    number of iterations of subproblem 
%          fX:    state of subproblem upon completion
%          rXVec: vector of residual outputs
%
%%
%% Code

             xFun  = @(x)vec( E.mtimes2(reshape(x,[Nx,r])*v')*v ) + lamx*x; % Optimising function            
[X,fX,rX,nX,rXvec] = minres(xFun,  vec(D*v + lamx*Xp), tol, Nin, [], [], vec(X));
             X     = reshape(X, Nx,r);
end