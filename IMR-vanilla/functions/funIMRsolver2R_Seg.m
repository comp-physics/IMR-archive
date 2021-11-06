function [t,R,U,tdel,Tdel,Cdel] = funIMRsolver2R_Seg(model,matprop,tspanStart,tspanEnd,R0,NT,NTM, ...
    Pext_type,Pext_Amp_Freq,disptime,Tgrad,Tmgrad,Cgrad,Dim,comp,Pfit_pp,Pfit_pp_d1,eqR)

% Authors:
% Carlos Barajas
% carlobar@umich.edu
% Umich Mechanical Engineering BS '16
% Jon Estrada
% jonathan_estrada@brown.edu
% Brown Solid Mechanics, PhD '17

% last Update: 8/21/2017

% Inputs:
% tspan - time to run simulation
% R0 - Initial Radii
% NT - number of nodes for temperature and concentration fields
% NTM - number of nodes for temperature in the medium
% Pext_type - type of external pressure ('sn' = sine, 'RC' = Rayleigh
% collapse, 'RG' = Rayleigh growth, impulse 'ip',...
% non-equlibrium initial conditions (laser caviation and Flynn(1975) ) 'IC'
% Pext_Amp_Freq - amplitude and frequency of external pressure [amp w]

% Note: For this code the out-of-equilibrium Rayleigh Collapse the intial
% mass in the bubble and radii are specified

% FOR THE FOLLOWING INPUTS 0 = FALSE AND 1 = TRUE
% disptime - Displays elapsed time on the command window
% Tgrad - Models temperature gradients in the bubble
% Tmgrad- Models temperature gradients outside the bubble
% Cgrad - Models concentration gradients in the bubble

% Outputs:
% t - time vector
% T_Bubble - Temperature inside the bubble
% T_Medium - Temperature outside the bubble
% R - Bubble Radius
% U - Bubble velocity
% P - Internal bubble pressure
% C - Vapor Concentration in the bubble
% Tm - Temperature in the medium
% Dim - outputs variables in dimensional form
% Comp - 0 (ignores compressibility effects) or 1 (uses Keller- Miksis)

%********************************************************************
% Citations for code:
%M.T. Warnez and E. Johnsen, "Numerical modeling of bubble dynamics in
%viscoelastic media with relaxation," Phys. Fluids 27, (2015).

%R. Gaudron, M.T. Warnez, and E. Johnsen, "Bubble dynamics in a
%viscoelastic medium with nonlinear elasticity,"
%J. Fluid Mech. 766, 54-75 (2015).

%X. Yang and C.C. Church, "A model for the dynamics of gas bubbles
%in soft tissue," J.Acoust. Soc. Am. 118, 3595-3606 (2005).

%A. Prosperetti, L. A. Crum, and K.W. Commander, "Nonlinear bubble
%dynamics," J.Acoust. Soc. Am. 83, 502-514 (1988).

%A.T. Preston, "Modeling heat and mass transfer in bubbly cavitating
%flows and shockwaves in cavitating nozzles," Ph.D. thesis,
%California Institute of Technology (2004).

%R.I. Nigmatulin, N.S. Khabeev, and F.B. Nagiev, "Dynamics, heat and mass
%transfer of vapour-gas plobubbles in a liquid," Int. J.
%Heat Mass Transfer, 24, 1033-1044 (1981).
%*************************************************************************

%*************************************************************************
% Assign values of material properties
G = matprop.G;
G1 = matprop.G1;
mu = matprop.mu;
alpha = matprop.alpha;
lambda_nu = matprop.lambda_nu;

%***************************************
% Load Parameters :
Pmt = IMRcall_parameters(R0,G,G1,mu); % Calls parameters script
k = Pmt(1); chi = Pmt(2); fom = Pmt(3); foh = Pmt(4); Ca = Pmt(5);
Re = Pmt(6); We = Pmt(7); Br = Pmt(8); A_star = Pmt(9); B_star = Pmt(10);
Rv_star = Pmt(11); Ra_star = Pmt(12); P0_star = Pmt(13); t0 = Pmt(14);
C0 = Pmt(15); L = Pmt(16); L_heat_star = Pmt(17); Km_star = Pmt(18);
P_inf = Pmt(19); T_inf = Pmt(20); C_star = Pmt(21); De = Pmt(22);
  
%****************************************

% Material Choice
fung = 0;
neoHook = 0;
nhzen = 0;
sls = 0;
linkv = 0;
if strcmp(model,'neoHook') == 1
    neoHook = 1;
elseif strcmp(model,'fung') == 1
    fung = 1;
elseif strcmp(model,'nhzen') == 1
    nhzen = 1;
elseif strcmp(model,'sls') == 1
    sls = 1;
elseif strcmp(model,'linkv') == 1
    linkv = 1;
else
    nhzen = 1;
end


%******************************************
% Create spatial nodes

% Inside the bubble
N = NT-1;
deltaY = 1/N;
i = 1:1:N+1;
yk = ((i-1)*deltaY)';

% Outside the bubble
Nm = NTM-1;
deltaYm = -2/Nm;
j = 1:1:Nm+1;
xk = (1+(j-1)*deltaYm)';
yk2 = ((2./(xk+1)-1)*L+1);

%******************************************
% Initial Conditions
tspanStart_star = tspanStart/t0; tspanEnd_star = tspanEnd/t0; 
R0_star = 1;
U0_star = 0;  % Change as needed
%Z10 = 0;
S0 = 0;


if strcmp(Pext_type,'ga')
    dt_star = Pext_Amp_Freq(2)/t0;
    tw_star = Pext_Amp_Freq(3)/t0;
end

% Need to modify intial conditions for the Out-of-Equilibrium Rayleigh
% Collpase:
% if strcmp(Pext_type,'IC')
%     Pv = Pvsat(1*T_inf)/P_inf; 
%     P0_star = Pext_Amp_Freq(1)/P_inf + Cgrad*Pvsat(1*T_inf)/P_inf;
%     % Need to recalculate intital concentration
%     
%     Pg000 = Pext_Amp_Freq(1)
%     % Calculate the equilibrium radii ratio for initial stress state:
%     [REq,~,~] = IMRCalc_Req(R0, Tgrad, Cgrad, Pext_Amp_Freq(1), G, G1, mu);
%     
%     REq
%      
%     %REq = 1; %removed 6/15/16 by Jon
%     C0 = C0*ones(1,NT);
%     %U0_star = -1*(1-P0_star)/(C_star); %Intitial velocity due to shockwave
%     U0_star = 0;
%     if sls == 1 || linkv == 1
%         S0 = -4/(3*Ca)*(1-REq^3);
%     elseif nhzen == 1 || neoHook == 1
%         S0 = -1/(2*Ca)*(5-REq^4-4*REq);
%     end
%     
% end
REq=eqR/R0;
S0 = -1/(2*Ca)*(5-REq^4-4*REq);

%X0 = [R0_star U0_star P0_star S0 Tau0 C0 Tm0];
X0 = [R0_star U0_star S0];

tau_del= [];
tdel=[];
Tdel = [];
Cdel = [];

%************************************************
% March equations in time
% options = odeset('RelTol',1e-10);
[t , X] = ode23tb(@bubble, [0, tspanEnd_star-tspanStart_star ] , X0);

R = X(:,1); % Bubble wall Radius
U = X(:,2); % Bubble wall velocity
S = X(:,3); % Stress integral
% P = X(:,3); % Internal pressure
% S = X(:,4); % Stress integral
% Tau = X(:,5:(NT+4)); % Variable relating to internal temp
% C =  X(:,(NT+5):(2*NT+4)); % Vapor concentration in the bubble
% Tm = X(:, (2*NT+5):end ); % Temperature variation in the medium
% T = (A_star -1 + sqrt(1+2*Tau*A_star)) / A_star; % Temp in bubble

% ******************************
% Transform variables back into their dimensional form
if (Dim == 1)
    R = R*R0;
    t = t*t0+tspanStart;
    % T = T*T_inf;
    % P = P*P_inf;
    S = S*P_inf;
    U = U*(R0/t0);
    tdel= tdel*t0;
    Tdel = Tdel*T_inf;
end
%***********************

%%
%*************************************************************************
% Nested function; ODE Solver calls to march governing equations in time
% This function has acess to all parameters above

    function dxdt = bubble(t,x)
        
        % Break x vector into indv. values
        R = x(1); % Bubble wall Radius
        U = x(2); % Bubble wall velocity
        S = x(3); % stress integral
%         P = x(3); % Internal pressure
%         %Z1 = x(4); % Stress integral
%         S = x(4);
%         Tau = x(5:(NT+4));
%         C = x((NT+5):(2*NT+4));
%         Tm = x((2*NT+5):end);
        
        % waitbar(t/(tspanEnd_star-tspanStart_star));
      
        % *********Solves for boundary condition at the wall**************
        %         if (Tmgrad == 1)
        %             if t/tspan_star> 0.001
        %                 %Might need to tune 0.001 for convergence:
        %                 guess= -.001+tau_del(end);
        %                 prelim  = fzero(@Boundary,guess);
        %             else
        %                 guess = -.0001;
        %                 prelim  = fzero(@Boundary,guess);
        %             end
        %         else
        prelim = 0;
        %        end
        
        %****************************************************************
%         % Sets value at boundary conditions
%         tau_del = [tau_del prelim];
%         Tau(end) = prelim;
%         T = TW(Tau);
%         Tm(1) = T(end);
%         % TW(prelim)
        
%         % Calculated variables
%         K_star = A_star*T+B_star;
%         C(end) =  CW(T(end),P);
%         
%         Rmix = C*Rv_star + (1-C)*Ra_star;
%         
%         % Gets variables that are not directly calculated as outputs
%         Tdel = [Tdel T(end)];
%         tdel = [tdel t];
%         Cdel = [Cdel C(end)];
        
        %Set external pressure
        %         if (Pext_type == 'sn')
        %             Pext =  -Pext_Amp_Freq(1)/P_inf*sin( Pext_Amp_Freq(2)*t*t0) ;
        %             P_ext_prime = -Pext_Amp_Freq(2)*t0*Pext_Amp_Freq(1)/P_inf...
        %                 *cos( Pext_Amp_Freq(2)*t*t0) ;
        %         elseif (Pext_type == 'RC')
        %             Pext = Pext_Amp_Freq(1)/P_inf ;
        %             P_ext_prime = 0;
        %         elseif (Pext_type == 'RG')
        %             Pext = -Pext_Amp_Freq(1)/P_inf ;
        %             P_ext_prime = 0;
        %         elseif (Pext_type == 'ip')
        %             Pext = -Pext_Amp_Freq(1)/P_inf*...
        %                 (1-heaviside(t-Pext_Amp_Freq(2)/t0)) ;
        %             P_ext_prime = 0;
        if strcmp(Pext_type,'IC')
            Pext = 0;
            P_ext_prime = 0;
        elseif strcmp(Pext_type , 'ga')
            Pext = pf(t)/P_inf;
            P_ext_prime = pfdot(t)/P_inf;
        end
        
        % *****************************************
        % Create derivative terms
        
%         % Temp. field of the gas inside the bubble
%         DTau  = D_Matrix_T_C*Tau;
%         DDTau = DD_Matrix_T_C*Tau;
%         
%         % Concentration of vapor inside the bubble
%         DC  = D_Matrix_T_C*C;
%         DDC = DD_Matrix_T_C*C;
%         
%         % Temp. field in the material
%         DTm = D_Matrix_Tm*Tm;
%         DDTm = DD_Matrix_Tm*Tm;
        
%         %***************************************
%         % Internal pressure equation
%         pdot = 3/R*(Tgrad*chi*(k-1)*DTau(end)/R - k*P*U +...
%             + Cgrad*k*P*fom*Rv_star*DC(end)/( T(end)*R* Rmix(end)* (1-C(end)) ) );
        % *****************************************
        
        %***************************************
%         % Temperature of the gas inside the bubble
%         U_vel = (chi/R*(k-1).*DTau-yk*R*pdot/3)/(k*P);
%         first_term = (DDTau.*chi./R^2+pdot).*(K_star.*T/P*(k-1)/k);
%         second_term = -DTau.*(U_vel-yk*U)./R;
%         
%         Tau_prime = first_term+second_term;
%         Tau_prime(end) = 0;
%         Tau_prime = Tau_prime*Tgrad;
        % *****************************************
        
        %***************************************
        % Vapor concentration equation
%         U_mix = U_vel + fom/R*((Rv_star - Ra_star)./Rmix).*DC;
%         one = DDC;
%         two = DC.*(DTau./(K_star.*T)+((Rv_star - Ra_star)./Rmix).*DC );
%         three =  (U_mix-U.*yk)/R.*DC;
%         
%         C_prime = fom/R^2*(one - two) - three;
%         C_prime(end) = 0;
%         C_prime = C_prime*Cgrad;
        %*****************************************
        
        %***************************************
        % Material temperature equations
%         first_term = (1+xk).^2./(L*R).*(U./yk2.^2.*(1-yk2.^3)/2+foh/R.*((xk+1)/(2*L)-1./yk2)).* DTm;
%         second_term = foh/R^2.*(xk+1).^4/L^2.*DDTm/4;
%         % third_term =  3*Br./yk2.^6.*(4/(3*Ca).*(1-1/R^3)+4.*U/(Re.*R)).*U./R;
%         % JY!!!
%         third_term = 3*Br./yk2.^6.*(4.*U/(Re.*R)).*U./R;
%         Tm_prime = first_term+second_term+third_term;
%         Tm_prime(end) = 0; % Sets boundary condition on temp
%         Tm_prime(1) = 0; % Previously calculated;
%         Tm_prime = Tm_prime*Tmgrad; %Tmgrad makes this quantity zero
        %*****************************************
        
         
        %***************************************
        % Elastic stress in the material
        if linkv == 1
            if  strcmp(Pext_type,'IC')
                Rst = R/REq;
                S = -4/(3*Ca)*(1 - 1/Rst^3) - 4/Re*U/R;
                Sdot =  -4/Ca*U/R/Rst^3 + 4/Re*U^2/R^2;
            else
                S = -4/(3*Ca)*(1 - 1/R^3) - 4/Re*U/R;
                Sdot =  -4/Ca*U/R^4 + 4/Re*U^2/R^2;
            end
        elseif neoHook == 1
            if  strcmp(Pext_type,'IC')
                Rst = R/REq;
                % ====== Old for neoHook ======
                S = -(5 - 4/Rst - 1/Rst^4)/(2*Ca) - 4/Re*U/R ;
                Sdot =  -2*U/R*(1/Rst + 1/Rst^4)/Ca + 4/Re*U^2/R^2; 
                % JY!!! "- 4/Re*udot/R" is added finally: Sdot = Sdot - SdotA*udot/R; % JY!!! Pay attention to here! 
            else
                S = -(5 -4/R - 1/R^4)/(2*Ca) - 4/Re*U/R;
                Sdot =  -2*U*(1/R^2 + 1/R^5)/Ca + 4/Re*U^2/R^2;
            end
        elseif fung == 1
            if strcmp(Pext_type,'IC')
                Rst = R/REq;
                % ====== JY!!! for Fung model ======
                % ******** JY!!! First order Fung model approx ======
                S = -(1-3*alpha)*(5 - 4/Rst - 1/Rst^4)/(2*Ca) - ...
                    2*alpha*(-27/40 - 1/8/Rst^8 - 1/5/Rst^5 -1/Rst^2 + 2*Rst)/(Ca) - 4/Re*U/R;
                Sdot = -2*U/R*(1-3*alpha)*(1/Rst + 1/Rst^4)/Ca - ...
                     2*alpha*U/R*(1/Rst^8 + 1/Rst^5 + 2/Rst^2 + 2*Rst)/(Ca) + 4/Re*U^2/R^2; 
                % JY!!! "- 4/Re*udot/R" is added finally: Sdot = Sdot - SdotA*udot/R; % JY!!! Pay attention to here! 
            else
                disp('not implemented.');
            end
                % ====== JY!!! First order Fung G + first order mu model approx ====== 
                % alpha = 0; lambda_nu = 0.001; 
%                 Lv = 1; 
%      
%                 zeta = linspace(-1,0.99,200);
%                  
%                 tempr = R*( (2./(1-zeta)-1)*Lv + 1 );
%                 tempr0 = (tempr.^3+R0^3-R^3).^(1/3);
%                 gammadot = -0.5*( 2*(tempr0.^2)./(tempr.^3) + 1./tempr0 ) *R^2./(tempr.^2) * U;
%                 % figure; plot(zeta,gammadot); pause;
%                 % tempmu = 1/Re .* heaviside(-abs(gammadot)+1/lambda_nu) .* (1-lambda_nu^2*(gammadot.^2));
%                 tempmu = 1/Re .* exp(-lambda_nu^2*(gammadot.^2));
%                 tempS = -12*2*tempmu*U/R .* (1-zeta).^2 ./ (2+(1-zeta)*(1/Lv-1)).^4 /(Lv^3);
%                  
%                 S = -(1-3*alpha)*(5 - 4/Rst - 1/Rst^4)/(2*Ca) - ...
%                     2*alpha*(-27/40 - 1/8/Rst^8 - 1/5/Rst^5 -1/Rst^2 + 2*Rst)/(Ca) + ...
%                     trapz(zeta,tempS);
%                  
%                 Sdot = -2*U/R*(1-3*alpha)*(1/Rst + 1/Rst^4)/Ca - ...
%                      2*alpha*U/R*(1/Rst^8 + 1/Rst^5 + 2/Rst^2 + 2*Rst)/(Ca); 
                
                % ******** JY!!! Second order Fung model approx ******
                % S = 2*(1-3*alpha+4.5*alpha^2)*(-5/4 + 1/Rst + 1/4/Rst^4)/(Ca) + ...
                %     2*(27/40*alpha-221/90*alpha^2 + alpha^2/24/Rst^12 + alpha^2/18/Rst^9 + (alpha-3*alpha^2)/8/Rst^8 + ... )
                %     2*alpha^2/6/Rst^6 + (alpha-3*alpha^2)/5/Rst^5 + 2*alpha^2/3/Rst^3 + (2*alpha-6*alpha^2)/2/Rst^2 - ...
                %     2*alpha^2*log(Rst) - (2*alpha-6*alpha^2)*Rst - 2/3*alpha^2*Rst^3)/(Ca) - ...
                %       4/Re*U/R;
                % Sdot = 2*U/R*(1-3*alpha+4.5*alpha^2)*(-1/Rst-1/Rst^4)/(Ca) + ...
                %     2*U/R*(-alpha^2/2/Rst^12 - alpha^2/2/Rst^9 - (alpha-3*alpha^2)/Rst^8 ...
                %     -2*alpha^2/Rst^6 - (alpha-3*alpha^2)/Rst^5 - 2*alpha^2/Rst^3 - (2*alpha-6*alpha^2)/Rst^2 ...
                %     -2*alpha^2 - (2*alpha-6*alpha^2)*Rst - 2*alpha^2*Rst^3)/(Ca) + 4/Re*U^2/R^2;

                % ******** JY!!! The following original Fung model codes don't work. ********
                % tempbeta = linspace(Rst,1,100);
                % tempS = -2*(tempbeta.^-5 + tempbeta.^-2) .* exp(alpha*(tempbeta.^-4+2*tempbeta.^2-3))/(Ca);
                % S = trapz(tempbeta,tempS) - 4/Re*U/R; 
                % S = sum( (1-Rst)/(length(tempbeta)-1) .* tempS );
                % Sdot = 2*U/REq*(Rst^5 + Rst^2)*exp(alpha*(Rst^(-4)+2*Rst^2-3))/(Ca) + ...
                %           4/Re*U^2/R^2;
        elseif sls == 1
            if  strcmp(Pext_type,'IC')
                Rst = R/REq;
                Sdot = -S/De - 4*(1-1/Rst^3)/(3*Ca*De) - 4/(Re*De)*U/R - 4*U/(Ca*R);
            else
                Sdot = -S/De - 4*(1-1/R^3)/(3*Ca*De) - 4/(Re*De)*U/R - 4*U/(Ca*R);
            end
        elseif nhzen == 1
            if  strcmp(Pext_type, 'IC')
                Rst = R/REq;
                Sdot = -S/De - 1/(2*Ca*De)*(5-1/Rst^4-4/Rst)-4*U/(R*Re*De)...
                    -4*U/(R*Ca)/(Rst^3-1)*(3/14*Rst^3+Rst^2-3/(2*Rst)+2/(7*Rst^4));
                if isinf(Sdot)
                    Rst=Rst+eps;
                    Sdot = -S/De - 1/(2*Ca*De)*(5-1/Rst^4-4/Rst)-4*U/(R*Re*De)...
                        -4*U/(R*Ca)/(Rst^3-1)*(3/14*Rst^3+Rst^2-3/(2*Rst)+2/(7*Rst^4));
                end
            else
                Sdot = -S/De - 1/(2*Ca*De)*(5-1/R^4-4/R)-4*U/(R*Re*De)...
                    -4*U/(R*Ca)/(R^3-1)*(3/14*R^3+R^2-3/(2*R)+2/(7*R^4));
                if isinf(Sdot)||isnan(Sdot)
                    R = R+eps;
                    Sdot = -S/De - 1/(2*Ca*De)*(5-1/R^4-4/R)-4*U/(R*Re*De)...
                        -4*U/(R*Ca)/(R^3-1)*(3/14*R^3+R^2-3/(2*R)+2/(7*R^4));
                end
            end
        end
        %****************************************************
        
        % Equations of motion
        rdot = U;
        
        %         if (Tgrad == 0)
        %             P = P0_star*(1/R)^(3*k);
        %             pdot = -3*k*U/R*P;
        %         end
        
        % JY!!! dangeraous!!!  Here is only true for Cgrad=1;
        Pv = 0;
        % Pv = (Pvsat(T(end)*T_inf)/P_inf);
        if comp == 0
            %Rayleigh-Plesset equation
            udot = (P + abs(1-Cgrad)*Pv  - 1 - Pext + S -1/(We*R) -1.5*U^2)/R;
        else
            % Keller-Miksis equation
            if linkv==1 || neoHook==1 || fung==1
                SdotA = 4/Re;
            elseif sls==1 || nhzen==1
                SdotA = 0;
            end
            % P_star = P/P_inf; Pdot_star = Pdot/P_inf*t0;
            P = ppval(Pfit_pp,t*t0+tspanStart)/P_inf; pdot = ppval(Pfit_pp_d1,t*t0+tspanStart)/P_inf*t0;
            
            if P>1, P =1; end; 
            if pdot>1 , pdot=1 ; end;
            
%             udot = Re*C_star/(Re*(C_star-U)*R+4)*( (1+U/C_star)*(P-1/We/R-1-1/2/Ca*(5-4*REq/R-(REq/R)^4)-4/Re*U/R) + ...
%                 R/C_star*(pdot+U/We/R^2-2*U/R/Ca*(REq/R+(REq/R)^4) + 4*U^2/Re/R^2 ) - ...
%                 1.5*(1-U/3/C_star)*U^2 );
            udot = ((1+U/C_star)...
                *(P  + abs(1-Cgrad)*Pv -1/(We*R) + S - 1 - Pext)  ...
                + R/C_star*(pdot+ U/(We*R^2) + Sdot -P_ext_prime ) ...
                - 1.5*(1-U/(3*C_star))*U^2)/((1-U/C_star)*R);%+JdotA/(C_star));
 
        end
        % ****************************************
        % ====== Keep viscosity the same ======
         Sdot = Sdot - SdotA*udot/R; % JY!!! Pay attention to here! 
        % ====== Change viscosity ======
%         temprdot = U * ( (2./(1-zeta) -1 )*Lv + 1 );
%         tempr0dot = (tempr.^2.*temprdot - R^2*U)./(tempr0.^2);
%         mudot_part0 = 2*(2*(tempr0.^2)./(tempr.^3) + 1./tempr0) * R^4 ./ tempr.^4 * U^2 .* ... 
%                       ( ( 4*tempr0.*tempr0dot.*tempr-6*tempr0.^2.*temprdot ) ./ (tempr.^4) - tempr0dot./(tempr0.^2) );
%         mudot_part2 = (2*tempr0.^2./tempr.^3 + 1./tempr0).^2 * U^2 .* (4*R^3*(U.*tempr-R.*temprdot))./(tempr.^5); 
%         mudot_part1 = (2*tempr0.^2./tempr.^3 + 1./tempr0).^2 *2*U * udot * R^4 ./ (tempr.^4); % JY!!!
%         % mudot = - (lambda_nu^2)/4/Re.*(mudot_part0+mudot_part1+mudot_part2) .* heaviside(-abs(gammadot)+1/lambda_nu);
%         mudot = exp(-(lambda_nu^2)*(gammadot.^2)) .* ( - (lambda_nu^2)/4/Re.*(mudot_part0+mudot_part1+mudot_part2) );
%          
%         tempSdot = -12*2*(1-zeta).^2 ./ ((2+(1-zeta)*(1/Lv-1)).^4*Lv^3) .* ...
%                    ( mudot*U*R + 1/Re*udot*R - 1/Re*U^2 ) / R^2;  
%          
%         Sdot = Sdot + trapz(zeta,tempSdot);
        % ====== End of this part ======
        % ****************************************
        
        % dxdt = [rdot; udot; pdot; Sdot; Tau_prime; C_prime; Tm_prime];
        dxdt = [rdot; udot; Sdot];
    end
%*************************************************************************

% Other nested functions used to carry out repetetive calculations
% throughout the code

%     function Tw = TW(Tauw)
%         %calculates the temperature at the bubble wall as a fuction of \tau
%         Tw = (A_star -1 + sqrt(1+2*Tauw*A_star)) / A_star;
%     end
% 
%     function Cw = CW(Tw,P)
%         % Calculates the concentration at the bubble wall
%         
%         %Function of P and temp at the wall
%         theta = Rv_star/Ra_star*(P./(Pvsat(Tw*T_inf)/P_inf) -1);
%         Cw = 1./(1+theta);
%     end
% 
%     function Tauw = Boundary(prelim)
%         % Solves temperature boundary conditions at the bubble wall
%         % Create finite diff. coeffs.
%         % Coefficients in terms of forward difference
%         
%         %    %Second order
%         coeff = [-3/2 , 2 ,-1/2 ];
%         Tm_trans = Tm(2:3);
%         T_trans = flipud(Tau(end-2:end-1));
%         C_trans = flipud(C(end-2:end-1));
%         
%         %   Could implement any order... sixth order example is shown below
%         
%         %     Sixth order
%         %     coeff= [-49/20 ,6	,-15/2	,20/3	,-15/4	,6/5	,-1/6]; %Sixth order coeff
%         %     Tm_trans = Tm(2:7);
%         %     T_trans = flipud(Tau(end-6:end-1));
%         %     C_trans = flipud(C(end-6:end-1));
%         
%         Tauw =chi*(2*Km_star/L*(coeff*[TW(prelim); Tm_trans] )/deltaYm) +...
%             chi*(-coeff*[prelim ;T_trans] )/deltaY + Cgrad*...
%             fom*L_heat_star*P*( (CW(TW(prelim),P)*(Rv_star-Ra_star)+Ra_star))^-1 *...
%             (TW(prelim) * (1-CW(TW(prelim),P))  ).^(-1).*...
%             (-coeff*[CW(TW(prelim),P); C_trans] )/deltaY;
%         %************************************************************************
%     end
% 
% % Gaussian pressure functions
% % acoustic pressure
%     function p = pf(t)
%         if t<(dt_star-5*tw_star) || t>(dt_star+5*tw_star)
%             p=0;
%         else
%             p = -Pext_Amp_Freq(1)*exp(-(t-dt_star).^2/tw_star^2);
%         end
%     end
% 
% % time derivative of acoustic pressure
%     function pdot = pfdot(t)
%         if t<(dt_star-5*tw_star) || t>(dt_star+5*tw_star)
%             pdot=0;
%         else
%             pdot = 2*(t-dt_star)/tw_star^2*Pext_Amp_Freq(1).*exp(-(t-dt_star).^2/tw_star^2);
%         end
%         
%     end

end




