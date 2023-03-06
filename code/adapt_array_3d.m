function [cmap]=adapt_array_3d(yn);

yn=permute(yn,[4,1,2,3]);
[nc,ny,nx,nz]=size(yn);
rn=eye(nc);

% Find coil with maximum intensity for phase correction
[mm,maxcoil]=max(sum(sum(sum(permute(abs(yn),[3 2 4 1])))));   

if nx~=nz
    bs1=8;  % x-block size
    bs2=8;  % y-block size
    bs3=4;  % z-block size
    st=4;   % Increase to set interpolation step size
    st_z=2;
else
    bs1=8;  
    bs2=8;  
    bs3=8;  
    st=4;   
    st_z=4;
end

cmapsmall=zeros(nc,round(ny./st),round(nx./st),round(nz./st_z));
for z=st_z:st_z:nz
for x=st:st:nx
for y=st:st:ny
    %Collect block for calculation of blockwise values
    ymin1=max([y-bs1./2 1]);                   
    xmin1=max([x-bs2./2 1]);  
    zmin1=max([z-bs3./2 1]); 
    % Cropping edges
    ymax1=min([y+bs1./2 ny]);                 
    xmax1=min([x+bs2./2 nx]); 
    zmax1=min([z+bs3./2 nz]);

    ly1=length(ymin1:ymax1);
    lx1=length(xmin1:xmax1);
    lz1=length(zmin1:zmax1);
    m1=reshape(yn(:,ymin1:ymax1,xmin1:xmax1,zmin1:zmax1),[nc,lx1*ly1*lz1]);
      
    m=m1*m1'; % Signal covariance
      
    % Eigenvector with max eigenvalue for optimal combination
    [e,v]=eig(inv(rn)*m);                    
                                               
    v=diag(v);
    [mv,ind]=max(v);
      
    normmf=e(:,ind);
    
    % Phase correction based on coil with max intensity
    normmf=normmf.*exp(-j*angle(normmf(maxcoil)));

    cmapsmall(:,y./st,x./st,z./st_z)=normmf;
end
end
end

% Interpolation of weights upto the full resolution
% Done separately for magnitude and phase in order to avoid 0 magnitude 
% pixels between +1 and -1 pixels.
for i=1:nc
    cmap(i,:,:,:)=resize(squeeze(abs(cmapsmall(i,:,:,:))),[ny,nx,nz],'*bilinear').*exp(j.*resize(squeeze(angle(cmapsmall(i,:,:,:))),[ny,nx,nz],'*nearest'));
end

cmap=permute(cmap,[2,3,4,1]);
