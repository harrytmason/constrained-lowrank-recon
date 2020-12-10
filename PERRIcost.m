function c = PERRIcost(d,E,u,v,lamu,lamv,Up,Vp)
%% The cost function of k-t PERRI
%
% Inputs:  d:   k-space data
%          E:   nufft structure + more 
%          u:   current guess of u
%          v:   current guess of v
%          lam: the weighting on the priors [lambda_u, lambda_v]
%          Vp:  the prior of u
%          Vp:  the prior of v
%
% Outputs: c:   the current estimate of the cost function
%
%%
%% Code

c =     norm(  vec(d - E*(u*v')), 'fro' )^2 + ... % data consistency (reshape only as it may end up in 3D during coils - gives same output)
   lamu^2*norm(  u - Up      , 'fro' )^2 + ... % u prior consistency
   lamv^2*norm(  v - Vp      , 'fro' )^2;      % v prior consistency
end