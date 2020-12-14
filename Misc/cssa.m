function [ssaVec ssaSum] = cssa(A,B,fW,T)
%Quick function which uses subspacea, but takes the cosine of the output
%(since that's nearly always what I want)
%
% fW allows for windowing to ensure no energy outside of possible radial
% sampling in k-space. Can be a boolean, or a window of appropriate size
%
% Check subspacea for A,B and T explanation
%
% Created by Harry Mason, 13/03/2019
% Check subspacea.m for actual details

% Check for inputs
if nargin<3; fW  = false; end

% Window Data
if numel(fW)>1; A = Window2Ddata(A,fW); % apply the window as an input (use if doing this function a lot)
                B = Window2Ddata(B,fW);
elseif fW;      A = Window2Ddata(A); 
                B = Window2Ddata(B); 
end

% Do cssa
if nargin<4, ssaVec = cos(subspacea(A,B));
else;        ssaVec = cos(subspacea(A,B,T));
end

% Check sum
ssaSum = sum(ssaVec);