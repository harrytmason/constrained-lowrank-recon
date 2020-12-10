function [u,rU,nU, fU,rUvec] = PERRIusub(D,E,u,v,lam,Up, tol, Nin, Nx,r)
%% Calcuates the solution to the u subproblem of k-t PERRI
%
% Inputs:  D:   nufft_adj(k-space data)
%          E:   nufft structure + more 
%          u:   current guess of u
%          v:   current guess of v
%          lam: lambda_u, the weighting on the prior
%          Up:  the prior of u
%          tol: tolerance of problem
%          Nin: number of iterations of subproblem  
%
% Outputs: u:     optimised u
%          rU:    residual output
%          nU:    number of iterations of subproblem 
%          fU:    state of subproblem upon completion
%          rUVec: vector of residual outputs
%
%%
%% Code
          %  sz    = size(u);
          %  cfU   = vec(D*v + lam*Up); % Cost function
             uFun  = @(x)vec( E.mtimes2(reshape(x,[Nx,r])*v')*v ) + lam*x; % Optimising function
            
[u,fU,rU,nU,rUvec] = minres(uFun,  vec(D*v + lam*Up), tol, Nin, [], [], vec(u));
             u     = reshape(u, Nx,r);
end