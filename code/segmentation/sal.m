%
% saliency map generation based on method described in:
%
%  "Salient Region Detection and Segmentation" by Achanta, et. al
%
% mjsottile@gmail.com // Dec. 2012
%
function m=sal(input_im)
    % dimensions
    [rows,cols,channels] = size(input_im);    

    if (channels == 3) then
        % step 1: go from sRGB color space to L*a*b*
        cform = makecform('srgb2lab');
        im = applycform(input_im,cform);
    else
        im = input_im;
    end
    
    % kernel size is determined by dimension of image with largest
    % extent.
    bigdim = max([rows,cols]);

    % bound the kernel sizes by 1/2 and 1/8th of the largest dimension
    minsize = bigdim/8;
    maxsize = bigdim/2;
    
    % range of kernel sizes
    range = maxsize-minsize;
    
    % create three kernel sizes, which span 25%, 50%, and 75% of the
    % size range.
    ksiz = minsize + [floor(range*0.25),floor(range*0.5),floor(range*0.75)];
    
    % determine amount of padding to put on the outside of the image.
    % this is to compensate for boundary issues for when the kernel
    % extends outside the image.
    hpad = ceil(max(ksiz) / 2);
    pad = hpad * 2;
    
    % create the new image as being full of pixels with the minimum
    % intensity from inside the original image
    imbig = ones(rows+pad,cols+pad) .* min(im(:));
    
    % put the original image in the middle of the large padded image
    imbig(hpad:(hpad+rows-1),(hpad:hpad+cols-1)) = im;
    
    % for each kernel size, create the convolution kernel and
    % compute the difference with the original image.
    % TODO: handle color channels properly.
    for i=1:length(ksiz)
        kern = ones(ksiz(i),ksiz(i)) / (ksiz(i)*ksiz(i));
        
        % note: for very large kernels, it is fastest to use an
        % FFT-based convolution implementation.  must use this instead
        % of conv2 from matlab for performance reasons.
        c{i} = abs(conv2fft(double(imbig),kern,'same') - double(imbig));
    end        
    
    % accumulate up the m vector
    mbig=zeros(rows+pad,cols+pad);
    for i=1:length(ksiz)
        mbig=mbig+c{i};
    end
    
    % extract the middle of the m matrix which corresponds to the
    % m-values for the original image
    m = mbig(hpad:(hpad+rows-1),(hpad:hpad+cols-1));
