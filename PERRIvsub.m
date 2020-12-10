function [v,rV,nV, fV,rVvec] = PERRIvsub(D,E,u,v,lam,Vp, tol, Nin,Nt,r)
%% Calcuates the solution to the v subproblem of k-t PERRI
%
% Inputs:  D:   nufft_adj(k-space data)
%          E:   nufft structure + more 
%          u:   current guess of u
%          v:   current guess of v
%          lam: lambda_v, the weighting on the prior
%          Vp:  the prior of v
%          tol: tolerance of problem
%          Nin: number of iterations of subproblem  
%
% Outputs: v:     optimised v
%          rV:    residual output of cost function
%          nV:    number of iterations of subproblem 
%          fV:    state of subproblem upon completion
%          rVVec: vector of residual outputs
%
%%
%% Code
             %sz    = size(v');
            % cfV   = vec(u'*D + lam*Vp');
             vFun  = @(x)vec(u'*E.mtimes2(u*reshape(x,[r, Nt]))) + lam*x;
             %vFun  = @(x)vec(reshape(x,[r, Nt])*reshape(x,[r, Nt])'*u'*E.mtimes2(u*reshape(x,[r, Nt]))) + lam*x;
[v,fV,rV,nV,rVvec] = minres(vFun, vec(u'*D + lam*Vp'), tol, Nin, [], [], vec(v'));
             v     = reshape(v, [r, Nt])';
end