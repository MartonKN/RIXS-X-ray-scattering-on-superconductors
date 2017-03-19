function plotH0(params)
    H0 = createH0(params);
    ix = 1:params.Nx; ix = ix - (params.Nx + 1)/2;
    iy = 1:params.Ny; iy = iy - (params.Ny + 1)/2;
    kx = (-pi):(2*pi/params.Nx):pi; kx = kx(1:params.Nx);
    ky = (-pi):(2*pi/params.Ny):pi; ky = ky(1:params.Nx);
    
    [IX,IY] = meshgrid(ix,iy);
    [KX,KY] = meshgrid(kx,ky);
    
    IX = reshape(IX,[numel(IX),1]);
    IY = reshape(IY,[numel(IY),1]);
    KX = reshape(KX,[numel(KX),1]);
    KY = reshape(KY,[numel(KY),1]);
    
    Ukr = zeros(length(IX),length(KX));
    for ir=1:length(IX)
        for ik=1:length(KX)
            Ukr(ir,ik) = exp(-1i*(IX(ir)*KX(ik) + IY(ir)*KY(ik)))/sqrt(params.Nx*params.Ny);
        end
    end
    
    KX = reshape(KX,[params.Nx,params.Ny]);
    KY = reshape(KY,[params.Nx,params.Ny]);
    
    Hk0 = diag(real(Ukr*H0*Ukr'));    
    Hk0 = reshape(Hk0,params.Nx,params.Ny);
    % surf(reshape(KX,[params.Nx,params.Ny]), reshape(KY,[params.Nx,params.Ny]), Hk0)
    figure(1); surf(KX, KY, Hk0)
    
    % As a cross check, let us determine the spectrum in Fourier space as
    % well
    Hk1 = zeros(length(kx), length(ky));
    [lJ, ~] = size(params.J);
    for ikx = 1:length(kx)
        for iky = 1:length(ky)
            for iJ = 1:lJ
                Dx = params.J(iJ,1);
                Dy = params.J(iJ,2);
                hopping = params.J(iJ,5);
                Hk1(ikx,iky) = Hk1(ikx,iky) - hopping * exp(1i*(kx(ikx)*Dx + ky(iky)*Dy));
            end
        end
    end
    Hk1 = real(Hk1);
    figure(2); surf(KX, KY, Hk1)
    figure(3); contour(KX, KY, Hk0,[0,0])
end