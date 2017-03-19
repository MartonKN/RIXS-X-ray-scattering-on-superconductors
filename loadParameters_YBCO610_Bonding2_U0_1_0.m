function params = loadParameters()
    Nx = 22; % #lattice sites in the x direction
    Ny = 22; % #lattice sites in the y direction
    Nb = 1;  % #bands
    neighbor_cutoff = 10; % cut off intensity after this number of sites
    
    T = 0.0070; % temperature (eV) 90 K = 0.0077 eV
    U0 = -1.000; % Core-hole potential (eV)
    Gamma = 0.250; % Inverse lifetime of the core hole (eV): 300-500 meV
    doping = 0.014; % Hole doping per Cu site
    apparatusResolution = 0.080; % resolution limit of the measurement
                                  % apparatus (eV). We convolve the final
                                  % signal with a Gaussian of this width.
    
    % Lattice symmetries
    SigmaXLatt = [1 , 0; 0, -1];
    InvLatt =    [-1, 0; 0, -1];
    SigmaDLatt = [0 , 1; 1,  0];
    symmetriesLatt = zeros(2,2,3);
    symmetriesLatt(:,:,1) = SigmaXLatt;
    symmetriesLatt(:,:,2) = InvLatt;
    symmetriesLatt(:,:,3) = SigmaDLatt;
    

    % List of all non-zero hopping amplitudes for Tl(2)Ba(2)CuO(6+delta)
    % t is a vector of the following lists: [Dx, Dy, mu, nu, tValue]

J00 = (-1)* 0.0500;
J10 = (-1)* -0.4200/4;
J11 = (-1)* 0.1163/4;
J20 = (-1)* -0.0983/4;
J21 = (-1)* 0.0353/8;
J22 = (-1)* 0.0/4;

    J = [];
    J = [J; [0, 0, 1,1, J00]];

    J = [J; [1, 0, 1,1, J10]];
    J = [J; [0, 1, 1,1, J10]];
    J = [J; [-1,0, 1,1, J10]];
    J = [J; [0,-1, 1,1, J10]];
    
    J = [J; [1, 1, 1,1, J11]];
    J = [J; [1,-1, 1,1, J11]];
    J = [J; [-1,1, 1,1, J11]];
    J = [J; [-1,-1,1,1, J11]];
    
    J = [J; [2, 0, 1,1, J20]];
    J = [J; [0, 2, 1,1, J20]];
    J = [J; [-2,0, 1,1, J20]];
    J = [J; [0,-2, 1,1, J20]];
    
    J = [J; [ 2, 1, 1,1, J21]];
    J = [J; [ 1, 2, 1,1, J21]];
    J = [J; [ 2,-1, 1,1, J21]];
    J = [J; [ 1,-2, 1,1, J21]];
    J = [J; [-2, 1, 1,1, J21]];
    J = [J; [-1, 2, 1,1, J21]];
    J = [J; [-2,-1, 1,1, J21]];
    J = [J; [-1,-2, 1,1, J21]];
    
    J = [J; [2, 2, 1,1, J22]];
    J = [J; [2,-2, 1,1, J22]];
    J = [J; [-2,2, 1,1, J22]];
    J = [J; [-2,-2,1,1, J22]];

    % Time resolution of the RIXS calculation
    % (David Benjamin set tFinal=30, dt=1.5, sFinal=45, ds=1.5. in his Mathematica file.)
    tFinal = 30; 
    dt     = 0.8;
    sFinal = 50;
    ds     = 1.0;
    
    % Frequency and momentum resolution
    omega = -2:(0.1/3):3;
    Domega = -0.2:0.02:2;
    
    qtmp = linspace(-pi,pi,Nx+1);
    qtmp = qtmp(1:Nx);
    [qX, qY] = meshgrid(qtmp, qtmp);
    q1 = [reshape(qX,[1,numel(qX)]); reshape(qY,[1,numel(qY)])];
     
    qtmp = linspace(0,2*pi,Nx+1);
    qtmp = qtmp(1:floor(Nx/2+1));
    q2 = [0*qtmp, qtmp; qtmp, qtmp];

    q =	[q1,q2];

    % Note: it is important that we choose wave vectors from the discrete
    % Brillouin zone, otherwise the simulation becomes non-reliable.
    
    params = struct('J', J, 'T', T, 'doping', doping, ...
                    'U0', U0, 'Gamma', Gamma,...
                    'apparatusResolution', apparatusResolution, ...
                    'Nx', Nx, 'Ny', Ny, 'Nb', Nb, ...
                    'symmetriesLatt', symmetriesLatt, ...
                    'neighbor_cutoff', neighbor_cutoff,...
                    'tFinal', tFinal, 'dt', dt, ...
                    'sFinal', sFinal, 'ds', ds, ...
                    'omega', omega, 'Domega', Domega, 'q', q);
end
