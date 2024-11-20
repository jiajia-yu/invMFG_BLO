function a = lsbd_dct(a)
%LSBD_DCT N-D 
%   last dimension single boundary condition discrete cosine transform.
%   B = LSBD_DCT(A) returns the discrete cosine transform of A 
%   with single homogeneous Neumann boundary condition on first dimension.
%   The array B is the same size as A and contains the
%   discrete cosine transform coefficients.
%
%   This function is modified from MIRT_DCT by 
%   Andriy Myronenko, see www.bme.ogi.edu/~myron
%   
%   See also LSBD_IDCT.

persistent siz ww ind isreala;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Check input
if (nargin == 0) || isempty(a) 
    error('Insufficient input');
end
isreala=isreal(a);
ndim=ndims(a);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Check for the row vector
transpose=0;
if (ndim==2) && (size(a,1)==1)
    transpose=1; a=a';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Check if the variable size has changed and we need to
%%% precompute weights and indicies
precompute=0;
if  ~exist('siz','var')
    precompute=1;
elseif abs(numel(siz)-ndims(a))>0
    precompute=1;
elseif sum(abs(siz-size(a)),2)>0
    precompute=1;
elseif isreala~=isreal(a)
    precompute=1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Precompute weights and indicies
if precompute
    siz=size(a);
    ndim=ndims(a);
    
    for i=1:ndim-1
        n=siz(i);
        
        ww{i} = 2*exp(((-1i*pi)/(2*n))*(0:n-1)')/sqrt(2*n);
        ww{i}(1) = ww{i}(1) / sqrt(2);
        ind{i}=bsxfun(@plus, [(1:2:n) fliplr(2:2:n)]', 0:n:n*(prod(siz)/n-1));
        if (siz(i)==1), break; end
    end
    
    i = ndim;
%     ww{i} = cos( pi/(4*siz(i)+2) * (1:2:2*siz(i)-1)' * (2*siz(i)-1:-2:1) );
    ww{i} = cos( pi/(4*siz(i)+2) * (1:2:2*siz(i)-1)' * (1:2:2*siz(i)-1) );
    ww{i} = ww{i} * (2/sqrt(2*siz(i)+1));
    ww{i} = ww{i}.';
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Actual multidimensional DCT. Handle 1D and 2D cases
%%% separately, because .' is much faster than shiftdim.
% check for 1D or 2D cases
if ndim==2
    if (min(siz)==1)
        a=trans(a,ww{1});if transpose, a=a'; end  % 1D case
    else
        a = trans(dct(a,ww{1},ind{1}).',ww{2}).'; % 2D case
%         a = dct(trans(a,ww{1}).',ww{2},ind{2}).'; % 2D case
    end
else
    % ND case (3D and higher)
    for i=1:ndim-1
        a=reshape(dct(reshape(a,siz(1),[]),ww{i},ind{i}), siz); % run dct along vectors (1st dimension)
        siz=[siz(2:end) siz(1)];                   % circular shift of size array by 1
        a=shiftdim(a,1);                           % circular shift dimensions by 1
    end
    i = ndim;
    a = reshape(trans(reshape(a,siz(1),[]),ww{i}),siz);
    siz = [siz(2:end) siz(1)];
    a = shiftdim(a,1);
end

end

function a = trans(a,ww)

a = ww*a;

end

function a=dct(a,ww,ind)
%DCT  Discrete cosine transform 1D (operates along first dimension)
isreala=isreal(a);
k=1; if ~isreala, ia = imag(a); a = real(a); k=2; end % check if complex

% k=1 (if 'a' is real) and 2 (if 'a' is complex)
for i=1:k
    a=a(ind);     % rearange
    a = fft(a);  % ifft
    a = real(bsxfun(@times,  ww, a)); % multiply weights
    % check if the data is not real
    if ~isreala
        if i==1, ra = a; a = ia; clear ia;  % proceed to imaginary part
        else, a = complex(ra,a); end         % finalize the output
    end
end
end


