function windowed_data = Window2Ddata(data,win)
%
% A function to take image-space data (with as many extra dimensions as you
% like), apply a 2d transform to convert to the fourier domain, window, and then jump back
%
% Input: data (x*y*anything you don't want windowed)
%        win (either a scalar radius, or a 2D window)
% Output: windowed_data (same dimensions as data)
if nargin <2
    win = 1; %Default assume a radius of 1 - although a 2d window is also allowed
end


Nd = size(data);

if numel(win)>1 %The user gave us a 2D window to work with - no need to generate one ourselves
    Nw = size(win);
    if Nd(1)==Nw(1)
        windowed_data = ifft2( ifftshift( fftshift( fft2( data ) ).*win ) );
        
    elseif sqrt(Nd(1)) == Nw(1)
        data = reshape(data,[sqrt(Nd(1)) sqrt(Nd(1)) Nd(2:end)]);
        windowed_data = ifft2( ifftshift( fftshift( fft2( data ) ).*win ) );
        windowed_data = reshape(windowed_data,Nd);
        
    else
        warning('Window size doesn''t match data size')
    end
    
else
    if Nd(1)==Nd(2) %preferred option
        Win = Window2D(Nd(1),win,0);
        windowed_data = ifft2( ifftshift( fftshift( fft2( data ) ).*Win ) );
        
    elseif sqrt(Nd(1)) == round(sqrt(Nd(1))) %If they've vectorised the first dimension
        data = reshape(data,[sqrt(Nd(1)) sqrt(Nd(1)) Nd(2:end)]);
        Win = Window2D(sqrt(Nd(1)),win,0);
        windowed_data = ifft2( ifftshift( fftshift( fft2( data ) ).*Win ) );
        windowed_data = reshape(windowed_data,Nd);
        
    else
        disp('Input 1''s first two dimensions aren''t equal, and dim 1 doesn''t have a square root. Are you sure the input data is the correct size?')
    end
end
end