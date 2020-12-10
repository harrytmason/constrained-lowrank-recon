function [ulamstr, vlamstr] = genLamStr(Urange,Vrange,spacing,zero_edge, inf_edge)
%% genLamStr.m
%
% Purpose: A simple function to create a 2 column matrix where each lambda
% pair is a unique combination within the xrange and yrange, and are
% ordered sensibly.
%
% Inputs:   Urange: a 2-element array with the 1e(xrange(1)) > 1e(xrange(2))
%           Vrange: as above, but the range of lambda_V rather than lambda_U
%           spacing: the spacing to use between log powers of U_range and V_range
%           zero edge: boolean of whether to put zeros around the lowest edge of U and V (if 2D, the first element is for u, the second is for v)
%           inf edge: boolean of whether to put infs around the highest edge of U and V (if 2D, the first element is for u, the second is for v)
%
% Outputs:  lambda. Col 1 = lambda_U values, Col 2 = lambda_V values
%
%%% Function by Harry Mason, 2018, University of Oxford

%% Initialisation
%
if nargin<5; inf_edge  = 0; end
if nargin<4; zero_edge = 0; end
if nargin<3; spacing   = 1; end

%% Creating lamU and lamV
%
if numel(Urange)>0
lamU = (Urange(1):spacing:Urange(end));
else
    lamU = [];
end
if numel(Vrange)>0
lamV = (Vrange(1):spacing:Vrange(end));
else
    lamV = [];
end

if isnan(sum(lamU)); lamU = []; end
if isnan(sum(lamV)); lamV = []; end

%% Final number of unique lambda_u and lambda_v values
%
Ulen = length(lamU);
Vlen = length(lamV);

ulamstr = cell(Ulen,1);
vlamstr = cell(Vlen,1);

%% Lambda values
for u = 1:Ulen
ulamstr{u} = strcat( '10^{' , num2str(lamU(u)) , '}' );

end

for v = 1:Vlen
vlamstr{v} = strcat( '10^{' , num2str(lamV(v)) , '}' );
end

%% Adding leading zeros?
%
if zero_edge(1)
    ulamstr = {'0' ulamstr{1:end}};
end
if zero_edge(end)
    vlamstr = {'0' vlamstr{1:end}};
end
if inf_edge(1)
    ulamstr = {ulamstr{1:end} '\infty'};
end
if inf_edge(end)
    vlamstr = {vlamstr{1:end} '\infty'};
end

%% End! Simple right?