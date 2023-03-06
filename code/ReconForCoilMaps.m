function ref=ReconForCoilMaps_XD(kdata,Traj3D,DensityCompen3D, Par,GPU);
Bas=size(kdata,1)/2;

if GPU
    [nx,nz,ntviews,nc]=size(kdata);
    Traj3D=reshape(Traj3D,[nx,nz*ntviews,3]);
    DensityCompen3D=reshape(DensityCompen3D,[nx,nz*ntviews]);
    kdata1=reshape(kdata,[nx,nz*ntviews,nc]);
    [nx,ntviews,nc]=size(kdata1);
    
    filter_factor=20;
    filter=kaiser(nx,filter_factor);
    kdata1=kdata1.*repmat(filter,[1,ntviews,nc]);
    kdata1=kdata1.*repmat(sqrt(DensityCompen3D),[1,1,nc]);
    kdata1=reshape(kdata1,[nx*ntviews,nc]);
    
    Img_Dim=[Bas,Bas,Par];
    param.E = SCGPUNUFFT2(reshape(Traj3D,[nx*ntviews,3]),reshape(DensityCompen3D,[nx*ntviews,1]),Img_Dim);
    clear ref;
    for ch=1:nc
        ref(:,:,:,ch)=param.E'*kdata1(:,ch);
    end
    ref=ref/max(abs(ref(:)));
else
    [nx,ntviews,nz,nc]=size(kdata);
    filter_factor=20;
    filter=kaiser(nx,filter_factor);
    kdata1=double(kdata.*repmat(filter,[1,ntviews,nz,nc]));
    kdata1=double(kdata1.*repmat(sqrt(DensityCompen3D),[1,1,nz,nc]));
    
    Img_Dim=[Bas,Bas,Par];
    param.E = MCNUFFT3D(double(Traj3D),double(DensityCompen3D),ones(Img_Dim));
    clear ref;
    for ch=1:nc
        ref(:,:,:,ch)=param.E'*kdata1(:,:,:,ch);
    end
    ref=ref/max(abs(ref(:)));
end
