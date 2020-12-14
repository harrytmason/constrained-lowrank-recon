function c = CLRcost(d,E,X,T,lamx,lamt,lams,Xp,Tp)
%% The cost function of k-t PERRI
%
% Inputs:  d:    k-space data
%          E:    nufft structure + more 
%          X:    current guess of X
%          T:    current guess of T
%          lamx: the weighting on the spatial constraint
%          lamt: the weighting on the temporal constraint
%          lams: the weighting on the smoothness constraint
%          Xp:   the prior of X
%          Tp:   the prior of T
%
% Outputs: c:   the current estimate of the cost function
%
%%
%% Code

c =       norm(  vec(d - E*(X*T')), 'fro' )^2 + ... % data consistency (reshape only as it may end up in 3D during coils - gives same output)
   lamx^2*norm(  X - Xp      ,      'fro' )^2 + ... % x prior consistency
   lamt^2*norm(  T - Tp      ,      'fro' )^2 + ... % t prior consistency
   lams^2*norm( diff(T,1,1)  ,      'fro' )^2;      % smoothness constraint
end