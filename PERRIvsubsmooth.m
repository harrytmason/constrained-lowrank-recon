function [v,rV,nV, fV,rVvec] = PERRIvsubsmooth(D,E,u,v,lamv,Vp, lams, tol, Nin,Nt,r)
%% Calcuates the solution to the v subproblem of k-t PERRI
%
% Inputs:  D:   nufft_adj(k-space data)
%          E:   nufft structure + more 
%          u:   current guess of u
%          v:   current guess of v
%          lam: [lambda_v lambda_s, the weighting on the prior and then the
%          smoothing strngth
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
             vFun  = @(x)vec(u'*E.mtimes2(u*reshape(x,[r, Nt]))) + lamv*x + lams * vec([zeros(r,1) diff(reshape(x,[r, Nt]),1,2) ]+[  fliplr(diff(fliplr(reshape(x,[r, Nt])),1,2)) zeros(r,1)]);
             %vFun  = @(x)vec(reshape(x,[r, Nt])*reshape(x,[r, Nt])'*u'*E.mtimes2(u*reshape(x,[r, Nt]))) + lam*x;
[v,fV,rV,nV,rVvec] = minres(vFun, vec(u'*D + lamv*Vp'), tol, Nin, [], [], vec(v'));
             v     = reshape(v, [r, Nt])';
end