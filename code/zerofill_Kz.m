function [kdata] = zerofill_Kz(kdata,z_dimension,partitions,slices,centerpartition,applyfilter)

currentsize=size(kdata,z_dimension);

% No zeropadding neeeded, so just return
if (currentsize==slices)
   return; 
end

% Determin the size of the zero-padded array
oldsize=size(kdata);
newsize=oldsize;
newsize(z_dimension)=slices;

% Preallocate memory for the zero-padded version
kdata_zf = zeros(newsize);

% Estimate the offset where the measure k-space data needs to be inserted 
% into the zero-padded array
% TODO: Verify that this is correct
shift=slices/2-centerpartition+1;

fprintf('Acquired partitions = %d, final slices = %d\n',currentsize, slices);

if applyfilter    
    % Calculate a filter that should be applied if partial fourier or 
    % reduced slice resolution has been used
    filterSlopeWidth=currentsize/8.;
    filter=ones(currentsize,1);
    for i = 0:filterSlopeWidth-1
        filter(i+1)          = sin( (i+1) * pi / (2. * filterSlopeWidth+1) );
        filter(currentsize-i)= sin( (i+1) * pi / (2. * filterSlopeWidth+1) );
    end

    fprintf('Applying partial Fourier / reduced slice resolution filter (slope %d)\n', filterSlopeWidth); 
    
    % Apply the filter
    for z = 1:currentsize
        kdata(:,:,z,:,:)=kdata(:,:,z,:,:).*filter(z);
    end
end

% Zero-pad the k-space data
kdata_zf(:,:,shift+1:shift+currentsize,:,:)=kdata;
kdata=kdata_zf;

% TODO: Find a way to make the insertion dimension depending on parameter
%       z_dimension, so that kdata with different order can be processed as
%       well.

% To display the projections:
% as( ifftshift( ifft( ifftshift( ifft( fftshift( fftshift( kdata(:,:,:,2),1),3) ,[],1) ,1) ,[],3 ) ,3))
