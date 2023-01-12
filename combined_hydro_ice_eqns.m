function F = combined_hydro_ice_eqns(QNShuxg,params)
    % unpack variables
    M = params.M;
    Q_in = params.Q_in;
    Nx = params.Nx;
    Nh = params.Nh;
    ice_start = params.ice_start;
    Q = QNShuxg(1:Nh);
    N = QNShuxg(Nh+1:2*Nh);
    S = QNShuxg(2*Nh+1:3*Nh);
    h = QNShuxg(ice_start+1:ice_start + Nx);
    u = QNShuxg(ice_start + Nx+1:ice_start + 2*Nx);
    xg = QNShuxg(ice_start + 2*Nx+1);
    hf = (-bed(xg.*params.x0,params)/params.h0)/(params.r);

    %grid params unpack
    dt = params.dt/params.t0;
    dth= params.dt/params.th0;
    ds = params.dsigma;
    Nx = params.Nx;
    N1 = params.N1;
    sigma = params.sigma;
    sigma_elem = params.sigma_elem;
    b = -bed(xg.*sigma.*params.x0,params)/params.h0;
    dbdx = gradient(bed(xg.*params.sigma_h.*params.x0,params))./gradient(params.sigma_h.*params.x0*xg);
    Fh = zeros(Nx,1);
    Fu = zeros(Nx,1);
    
    %physical params unpack
    m     = 1/params.n;
    nglen = params.n;
    accum = params.accum;
    a = accum/params.a0;
    gamma = params.gamma;
    alpha = params.alpha;
    ss = params.transient;
       
    %previous time step unpack
    h_old = params.h_old;
    xg_old = params.xg_old;
    
    % Assign values based off coupling parameters
    if params.hydro_u_from_ice_u
        u_ice_interp = interp1(sigma,u,params.sigma_h,"linear","extrap");
    else
        u_ice_interp = 1000.*ones(params.Nx,1)./(params.year*params.u0);
    end

    if params.hydro_psi_from_ice_h
        %gradient_h = zeros(size(sigma_elem));
        %gradient_h(1:N1) = gradient(h(1:N1).*params.h0)./gradient(sigma_elem(1:N1).*params.x0*xg);
        %gradient_h(N1+1:end) = gradient(h(N1+1:end).*params.h0)./gradient(sigma_elem(N1+1:end).*params.x0*xg);
        %gradient_h_interp = interp1(sigma_elem,gradient_h,params.sigma_h,'linear','extrap');
        %params.psi = (-params.rho_w*params.g.*dbdx-params.rho_i*params.g.*gradient_h_interp)./params.psi0;
        h_interp = interp1(sigma_elem,h.*params.h0,params.sigma_h,'linear','extrap');
        params.psi = (-params.rho_w*params.g.*dbdx-params.rho_i*params.g.*gradient(h_interp)./gradient(params.sigma_h.*params.x0*xg))./params.psi0;
    else
        params.psi = 1*(1-3*exp(-20.*params.sigma_h));
    end

    if params.ice_N_from_hydro
        N_ice = interp1(params.sigma_h,N,sigma);
    else
        N_ice = 100000*ones(size(h_old))./params.N0;
    end

    % Q
    fq(1) = Q(1) - Q_in; % Boundary condition

    fq(2:Nh-1) = (params.eps_r*(params.r-1))*abs(Q(2:Nh-1)).^3./S(2:Nh-1).^(8/3) +...
                        params.eps_r.*S(2:Nh-1).*N(2:Nh-1).^3 + M + ...
                        params.eps_r.*params.beta.*u_ice_interp(2:Nh-1).*(S(3:Nh)-S(1:Nh-2))./(2*xg*params.dsigma_h(2:Nh-1)) - ...
                        (Q(3:Nh)-Q(1:Nh-2))./(2*xg*params.dsigma_h(2:Nh-1));
    fq(Nh) = (params.eps_r*(params.r-1))*abs(Q(Nh)).^3./S(Nh).^(8/3) +...
                        params.eps_r*S(Nh).*N(Nh).^3 + M + ...
                        params.eps_r.*params.beta.*u_ice_interp(Nh).*(S(Nh)-S(Nh-1))./(xg*params.dsigma_h(Nh-1)) - ...
                        (Q(Nh)-Q(Nh-1))./(xg*params.dsigma_h(Nh-1)); % one sided difference instead
    % N 
    fn(1) = Q(1).*abs(Q(1))./(S(1).^(8/3)) - params.psi(1) - ...
                        params.delta*(N(2)-N(1))./(xg*params.dsigma_h(1)); % use 1 sided difference instead of symmetry argument
    fn(2:Nh-1) =  Q(2:Nh-1).*abs(Q(2:Nh-1))./(S(2:Nh-1).^(8/3)) - params.psi(2:Nh-1) - ...
                        params.delta*(N(3:Nh)-N(1:Nh-2))./(2*xg*params.dsigma_h(2:Nh-1));
    fn(Nh) = N(Nh) - params.N_terminus; % Boundary condition
    % S
    fs(1) = abs(Q(1)).^3./(S(1).^(8/3)) - ... 
                        S(1).*N(1).^3 + ...
                        (ss.*params.sigma_h(1)*(xg-xg_old)/dth - params.beta.*u_ice_interp(1)).*(S(2)-S(1))./(xg*params.dsigma_h(1)) - ...
                        ss.*(S(1)-params.S_old(1))./dth; % one sided difference
    fs(2:Nh-1)= abs(Q(2:Nh-1)).^3./(S(2:Nh-1).^(8/3)) - ... 
                        S(2:Nh-1).*N(2:Nh-1).^3 + ...
                        (ss.*params.sigma_h(2:Nh-1).*(xg-xg_old)./dth - params.beta.*u_ice_interp(2:Nh-1)).*(S(3:Nh)-S(1:Nh-2))./(2*xg*params.dsigma_h(2:Nh-1)) - ...
                        ss.*(S(2:Nh-1)- params.S_old(2:Nh-1))./dth; 
    fs(Nh)= abs(Q(Nh)).^3./(S(Nh).^(8/3)) - ... 
                        S(Nh).*N(Nh).^3 + ...
                        (ss.*params.sigma_h(Nh).*(xg-xg_old)./dth - params.beta.*u_ice_interp(Nh)).*(S(Nh)-S(Nh-1))./(xg*params.dsigma_h(Nh-1)) - ...
                        ss.*(S(Nh)-params.S_old(Nh))./dth; % one sided difference

    % ice sheet equations
    Fh(1)      = ss.*(h(1)-h_old(1))./dt + (2.*h(1).*u(1))./(ds(1).*xg) - a;
    Fh(2)      = ss.*(h(2)-h_old(2))./dt -...
                    ss.*sigma_elem(2).*(xg-xg_old).*(h(3)-h(1))./(2*dt.*ds(2).*xg) +...
                        (h(2).*(u(2)+u(1)))./(2*xg.*ds(2)) -...
                            a;
    Fh(3:Nx-1) = ss.*(h(3:Nx-1)-h_old(3:Nx-1))./dt -...
                    ss.*sigma_elem(3:Nx-1).*(xg-xg_old).*(h(4:Nx)-h(2:Nx-2))./(2*dt.*ds(3:Nx-1).*xg) +...
                        (h(3:Nx-1).*(u(3:Nx-1)+u(2:Nx-2)) - h(2:Nx-2).*(u(2:Nx-2)+u(1:Nx-3)))./(2*xg.*ds(3:Nx-1)) -...
                            a;
    
    Fh(N1) = (1+0.5*(1+(ds(N1)/ds(N1-1))))*h(N1) - 0.5*(1+(ds(N1)/ds(N1-1)))*h(N1-1) - h(N1+1);
                        
    Fh(Nx)     = ss.*(h(Nx)-h_old(Nx))./dt -...
                    ss.*sigma_elem(Nx).*(xg-xg_old).*(h(Nx)-h(Nx-1))./(dt.*ds(Nx-1).*xg) +...
                        (h(Nx).*(u(Nx)+u(Nx-1)) - h(Nx-1).*(u(Nx-1)+u(Nx-2)))./(2*xg.*ds(Nx-1)) -...
                            a;

	%velocity
    Fu(1)      = (alpha).*(1./(xg.*ds(1)).^((1/nglen)+1)).*...
                 (h(2).*(u(2)-u(1)).*abs(u(2)-u(1)).^((1/nglen)-1) -...
                  h(1).*(2*u(1)).*abs(2*u(1)).^((1/nglen)-1)) -...
                  gamma.*params.shear_scale.*N_ice(1).* (u(1)./(u(1) + N_ice(1).^nglen)).^m -...
                  0.5.*(h(1)+h(2)).*(h(2)-b(2)-h(1)+b(1))./(xg.*ds(1));
    Fu(2:Nx-1) = (alpha).*(1./(xg.*ds(2:Nx-1)).^((1/nglen)+1)).*...
                 (h(3:Nx).*(u(3:Nx)-u(2:Nx-1)).*abs(u(3:Nx)-u(2:Nx-1)).^((1/nglen)-1) -...
                  h(2:Nx-1).*(u(2:Nx-1)-u(1:Nx-2)).*abs(u(2:Nx-1)-u(1:Nx-2)).^((1/nglen)-1)) -...
                  gamma.*params.shear_scale.*N_ice(2:Nx-1).*(u(2:Nx-1)./(u(2:Nx-1) + N_ice(2:Nx-1).^nglen)).^m -...
                  0.5.*(h(2:Nx-1)+h(3:Nx)).*(h(3:Nx)-b(3:Nx)-h(2:Nx-1)+b(2:Nx-1))./(xg.*ds(2:Nx-1));
    Fu(N1)     = (u(N1+1)-u(N1))/ds(N1) - (u(N1)-u(N1-1))/ds(N1-1);
    Fu(Nx)     = alpha.*(1./(xg.*ds(Nx-1)).^(1/nglen)).*...
                 (abs(u(Nx)-u(Nx-1)).^((1/nglen)-1)).*(u(Nx)-u(Nx-1)) - 0.5*(1-params.r)*hf;
             
    Fxg        = 3*h(Nx) - h(Nx-1) - 2*hf; 

    F = [fq';fn';fs';Fh;Fu;Fxg];
