%% Numerical simulation of femtosecond pulse propagation in the presence of
% third-order nonlinearity and pulse compression by self-phase modulation
% Leonardo Redaelli

format long

% Parameters for which we want to solve the NL Schroedinger equation:
% Initial condition: NO CHIRP, CEP=0 (my choice, I will test different values)
lambda0 = 0.8; % [microm] central wavelength of pulse (carrier)
dt0 = 20; % [fs] temporal duration (FWHM) of pulse (envelope)
Ipeak = 1e12; % [W/cm^2] peak intensity Ipeak = 2/Z0*n0*|E|^2 since convention Electric Field(z,t)=E(z,t)*e^i(w0t-k0z)+c.c.
CEP = 0; % [rad] initial carrier envelope phase (I have notice a linear sensitivity of CEPout to CEPin in all scenarios)
L = 1; % [m] propagation length of the pulse in the nonlinear medium (Chi3)
N = 1e3; % number of steps of the Split-Step method (propagation direction mesh)
dl = L/N; % [m] step length (propagation direction mesh)
Tmax = 200; % [fs] maximum time, -Tmax/2 to +Tmax/2 (electric field mesh)
M = 1e4; % number of sampling points (electric field mesh)


% Saving the coefficients for Sellmeier equation in the form:
% n^2(lambda) = A + sum(B*lambda^2/(lambda^2+C))
% as [A,B1,C1,B2,C2,...], where
% [A] = dimensionless, [B] = microm^-2, [C] = microm^2
% all valid (at the least) from 0.5 microm up to 1.5 microm

% Dictionaries are available from Matlab 2022b. Despite that, it couldn't
% be done as I thought with dictionaries because they are different from
% Python dictionaries, they allow the "value" to be only a scalar
% n0 = dictionary('Air', coeff_Air, 'Ar', coeff_Ar, 'Fused Silica', ...);
% ERROR: coeff_Air = [n0_Air, n2_Air] is not a scalar

% https://refractiveindex.info/?shelf=other&book=air&page=Ciddor
n0_Air = [1, 14926.44e-8, 19.36e-6, 41807.57e-8, 7.434e-3];
% n0_Air = @(x) (1+14926.44e-8./(1-19.36e-6./x.^2)+41807.57e-8./(1-7.434e-3./x.^2)).^0.5;

% https://refractiveindex.info/?shelf=main&book=Ar&page=Borzsonyi
n0_Ar = [1, 20332.29e-8, 206.12e-6, 34458.31e-8, 8.066e-3];
% n0_Ar = @(x) (1+20332.29e-8./(1-206.12e-6./x.^2)+34458.31e-8./(1-8.066e-3./x.^2)).^0.5;

% https://refractiveindex.info/?shelf=glass&book=fused_silica&page=Malitson
n0_Fused_Silica = [1, 0.6961663, 0.0684043^2, 0.4079426, 0.1162414^2, 0.8974794, 9.896161^2];
% n0_Fused_Silica = @(x) (1+0.6961663./(1-(0.0684043./x).^2)+0.4079426./(1-(0.1162414./x).^2)+0.8974794./(1-(9.896161./x).^2)).^0.5;


% Saving the nonlinear refractive index coefficients
% https://opg.optica.org/oe/fulltext.cfm?uri=oe-18-25-25847&id=208397
n2_Air = 5.7e-23; % [n2] = m^2/W @ 800 nm and 1 atm
n2_Ar = 1.94e-22; % [n2] = m^2/W @ 800 nm and 1 atm
% https://www.rp-photonics.com/nonlinear_index.html
% n2_Air = 1.22e-22; % [n2] = m^2/W @ 308 nm and 1 atm
n2_Fused_Silica = 2.19e-20; % [n2] = m^2/W @ 1030 nm


% % DEBUG & FORMULAE TEST (DIMENSIONALITY CHECK)
% n0(n0_Air, 0.8) % same as website (refractiveindex.info)
% c0 = 3e8; % [m/s]
% -d2n0(n0_Air, 0.8)*0.8/c0 * 10^12 % [s/(microm*m)] and I want it in [ps/(km*nm)]
% % to compare with website => 10^12. I get -0.063149 website -0.062493
% -d2n0(n0_Air, 1.3)*1.3/c0 * 10^12 % I get -0.014365 website -0.012819
% d2k0(n0_Air, 0.8) * 10^3 % [fs^2/microm] and I want it in [fs^2/mm] to
% % compare with website => 10^3. I get 0.021471 website 0.021233.
% % Quite different. I didn't use any approximation. Might be an error of
% % mine since I have not been much aware of floating point errors during all
% % the summations that take place, or it might be the website that
% % approximated the formula. Nevertheless, it opens the possibility to use 
% % the program for many different wavelengths, so I'll continue with this.


n = 0;
while n == 0
    NLM = input("Select the nonlinear medium you want among the following " + ...
    "list by typing its correspondent number (e.g. 1): \n 1) Air, " + ...
    "\n 2) Ar, \n 3) Fused Silica\n");

    switch NLM
        case 1
            n = [n0_Air, n2_Air];
        case 2
            n = [n0_Ar, n2_Ar];
        case 3
            n = [n0_Fused_Silica, n2_Fused_Silica];
        otherwise
            disp('Inexistent medium selected! Try to press 1 and then press Enter')
    end
end


Z0 = 377; % [Ohm] vacuum impedance
v = 299792458; % [m/s] speed of light in vacuum
E0 = sqrt(Ipeak*1e4*Z0*0.5/n0(n(1:end-1), lambda0)); % [V/m] E(z=0,t=0) "electric field"
w0 = 2*pi*v/(lambda0*1e-6) * 1e-15; % [rad/fs] pulsation of the electric field (carrier)

% playing around to generate time and frequency axis
t = Tmax * ((0:1/M:1)-0.5); % [fs] time axis
f0 = w0/(2*pi);
fmax = (M+1)/Tmax; % [fs^-1] maximum sampling frequency (electric field mesh)
df = 1/Tmax; % [fs^-1] frequency sampling step (electric field mesh)
if mod(length(t), 2) == 1 % if odd number of points in the mesh
    half_df = df/2;
    f = -fmax/2+half_df:df:fmax/2-half_df; % [fs^-1] frequency axis
else % if even number of points in the mesh
    f = -fmax/2:df:fmax/2-df; % [fs^-1] frequency axis
end

sigmat0 = dt0/sqrt(8*log(2)); % [fs] sigmat0^2=FWHM^2/(8*ln2)
tau = t./sigmat0; % [dimensionless] time axis normalized to the sigmat0
f_norm = f.*sigmat0; % [dimensionless] frequency axis normalized to the sigmat0
f0_norm = f0.*sigmat0; % [dimensionless] carrier frequency normalized to the sigmat0


E_field = E0 * exp(-0.25*tau.^2) .* exp(1i*(w0*t+CEP)); % E(z,tau) "electric field" (w0*t)=([w0*sigmat0]*[t/sigmat0])=(w0_norm*tau)

U = exp(-0.25*tau.^2) .* exp(1i*(w0*t+CEP)); % U=E(z,tau)/|E0| normalized "electric field"
%U = sech(tau.^2).^2 .* exp(1i*(w0*t+CEP));
U_prop = zeros(N, length(U)); % save in memory U for all the propagation steps (Split-Step Method)
U_prop(1,:) = U;

%% Plots
% plot the initial Intensity (instantaneous and of the envelope)
figure(1)
plot(t, 0.5/Z0*n0(n(1:end-1), lambda0)*(E_field+conj(E_field)).^2,'LineWidth',2)
hold on
plot(t, 0.5/Z0*n0(n(1:end-1), lambda0)* (E_field.*conj(E_field)),'LineWidth',2)
title('Intensity before entering the nonlinear medium')
xlabel('Time in the pulse reference frame $t[fs]$','interpreter','latex') 
ylabel('Intensity $I[W/m^2]$','interpreter','latex')
legend('Instantaneous Intensity','Averaged Intensity')
hold off

% plot the Initial Electric Field (Electric Field(z,t)=E(z,t)+c.c.)
figure(2)
plot(t, E_field+conj(E_field),'LineWidth',2)
title('Electric Field E before entering the nonlinear medium')
xlabel('Time $t[fs]$','interpreter','latex') 
ylabel('Electric Field $E[V/m]$','interpreter','latex')


%% Applying the Split-Step Method!
% Expressing the NL Schroedinger equation in the Wikipedia Split-Step
% method notation to solve it: https://en.wikipedia.org/wiki/Split-step_method
beta2 = -sign(d2k0(n(1:end-1), lambda0))/Ld(n(1:end-1), lambda0, dt0);
g = -1/Ln(n0(n(1:end-1), lambda0), n(end), lambda0, E0);
disp("Dispersion Length: " + num2str(Ld(n(1:end-1), lambda0, dt0)) + "[m]")
disp("Nonlinear Length: " + num2str(Ln(n0(n(1:end-1), lambda0), n(end), lambda0, E0)) + "[m]")

% The following works thanks to the fact that the NL schroedinger equation
% have an exact solution in the two separate cases:
% - when placing to 0 the dispersion term -> temporal domain
% - when placing to 0 the nonlinear term -> frequency domain
% Otherwise one would have needed to solve the eq. with a solver to make
% the step

% %% 1ST ORDER SOLUTION -> L+N
% for k = 1:(N-1)
%     U_n = exp(1i*g*U_prop(k,:).*conj(U_prop(k,:))*dl) .* U_prop(k,:); % Nonlinear Step
%     F_U = exp(1i*beta2*0.5*(2*pi.*(f_norm-f0_norm)).^2 *dl) .* fftshift(fft(U_n)); % Linear Step
%     U_prop(k+1,:) = ifft(ifftshift(F_U));
% end
%% 2ND ORDER SOLUTION -> L/2+N+L/2
% from what I have verified the accuracy of this method is nearly the same
% with respect to the previous one
for k = 1:(N-1)
    F_U1 = exp(1i*beta2*0.25*(2*pi.*(f_norm-f0_norm)).^2 *dl) .* fftshift(fft(U_prop(k,:))); % First Half Linear Step
    U_n = exp(1i*g*(U_prop(k,:).*conj(U_prop(k,:)))*dl) .* ifft(ifftshift(F_U1)); % Nonlinear Step
    F_U2 = exp(1i*beta2*0.25*(2*pi.*(f_norm-f0_norm)).^2 *dl) .* fftshift(fft(U_n)); % Second Half Linear Step
    U_prop(k+1,:) = ifft(ifftshift(F_U2));
end



%% Convert back to Electric Field, calculate Spectrum, Spectral Phase and Dispersion Compensated Pulse
E_prop = U_prop*E0; % E(z,tau) "electric field" where tau=t-z/vg or (t-z/vg)/sigmat0, changes only in the plot
S_w_prop = fftshift(fft(ifftshift(E_prop,2), [], 2),2); % S(z,w) = |Spectrum| * e^i(Spectral Phase)

S_phase = unwrap(angle(S_w_prop(end,:)));
% - sign necessary to be coherent with the fact that material with
% GDD>0 (normal disperison) introduces effectively a GDD>0???????????????
% NO! it's correct as it is, without the - sign. This is due to the fact
% that I represented the pulse as e^[i(w0*t-k0*z)], and therefore it's
% Fourier Transform is e^[-iz(k0+k0'*(w-w0)+k0''*(w-w0)^2+...)], so to find
% the GDD one only has to take the opposite signs of the coefficients of
% the fitted polynomial!!!!!!
E_phase = unwrap(angle(E_prop(end,:)));

% Fitting the Spectral Phase with a polynomial up to the P^th order
P = 5;
coeff_fit = -PhaseFit(2*pi*f', 2*pi*f0, P, S_phase', S_w_prop(end,:));
S_phase_fit = BuildPhaseFit(2*pi*f', 2*pi*f0, -coeff_fit);


% Remember that GDD, TOD, etc.. are defined as [d^n(phi(w))/dw^n]@w0
disp("GDD: " + num2str(coeff_fit(3)*factorial(2)) + "[fs^2/rad]")
disp("TOD: " + num2str(coeff_fit(4)*factorial(3)) + "[fs^3/rad^2]")
disp("4^TH ORDER DISPERSION: " + num2str(coeff_fit(5)*factorial(4)) + "[fs^4/rad^3]")

% Compensate dispersion of the pulse by applying an opposite spectral phase
% with respect to the fitted one (using the fitted one is as if we measured it)
E_ftl = DispersionCompensation(S_w_prop(end,:), -S_phase);


%% Plots
figure(3)
yyaxis left
plot(t, 0.5/Z0*n0(n(1:end-1), lambda0)* (E_prop(end,:).*conj(E_prop(end,:))),'LineWidth',2)
hold on
plot(t, 0.5/Z0*n0(n(1:end-1), lambda0)* (E_prop(1,:).*conj(E_prop(1,:))),'LineWidth',2)
plot(t, 0.5/Z0*n0(n(1:end-1), lambda0)* (E_ftl.*conj(E_ftl)),'LineWidth',2)
title('Temporal Intensity and Phase after propagation')
xlabel('Time $t[fs]$','interpreter','latex')
ylabel('Temporal Intensity $I[W/m^2]$','interpreter','latex')
hold off
yyaxis right
plot(t, E_phase,'LineWidth',2)
ylabel('Phase $\Phi[rad]$','interpreter','latex')
legend('Pulse Intensity after propagation', 'Pulse Intensity before propagation', 'Pulse Intensity after prop. FTL', 'Phase')
%legend('Pulse Intensity after propagation', 'Pulse Intensity before propagation', 'Phase')

figure(4)
imagesc(t, 0:dl:L, E_prop+conj(E_prop))
colorbar
title('Electric Field E (pulse reference frame)')
xlabel('Time $t[fs]$ in the pulse reference frame','interpreter','latex') 
ylabel('Propagation distance $z[m]$','interpreter','latex')

% % imagesc(t, 0:dl:L, 0.5/Z0*n0(n(1:end-1), lambda0)*(E_prop.*conj(E_prop)))
% imagesc(t/1e3, (0:dl:L)/1e3, 0.5/Z0*n0(n(1:end-1), lambda0)*(E_prop.*conj(E_prop)))
% colorbar
% title('Intensity I (pulse reference frame)')
% % xlabel('Time $t[fs]$ in the pulse reference frame','interpreter','latex') 
% % ylabel('Propagation distance $z[m]$','interpreter','latex')
% xlabel('Time $t[ps]$ in the pulse reference frame','interpreter','latex') 
% ylabel('Propagation distance $z[km]$','interpreter','latex')

figure(5)
yyaxis left
plot(f, 0.5/Z0*n0(n(1:end-1), lambda0)* (S_w_prop(end,:).*conj(S_w_prop(end,:))),'LineWidth',2)
hold on
plot(f, 0.5/Z0*n0(n(1:end-1), lambda0)* (S_w_prop(1,:).*conj(S_w_prop(1,:))),'LineWidth',2)
hold off
title('Spectral Intensity and Spectral Phase after propagation')
xlabel('Frequency $f[fs^{-1}=PHz]$','interpreter','latex')
ylabel('Spectral Intensity $S[a.u.]$','interpreter','latex')
yyaxis right
plot(f, S_phase, f, S_phase_fit,'LineWidth',2)
ylabel('Phase $\Phi[rad]$','interpreter','latex')
legend('Final Spectral Intensity', 'Initial Spectral Intensity', 'Spectral Phase (final)', 'Fitted Spectral Phase (final)')

figure(6)
imagesc(f, 0:dl:L, sqrt(S_w_prop.*conj(S_w_prop)))
colorbar
title('Spectrum |S(w,z)|')
xlabel('Frequency $f[fs^{-1}=PHz]$','interpreter','latex') 
ylabel('Propagation distance $z[m]$','interpreter','latex')



% Energy conservation check:
Energy0 = sum(0.5/Z0*n0(n(1:end-1), lambda0)*(U_prop(1,:)+conj(U_prop(1,:))).^2); % [J/m^2]
Energy1 = sum(0.5/Z0*n0(n(1:end-1), lambda0)*(U_prop(end,:)+conj(U_prop(end,:))).^2); % [J/m^2]


%% Save figures
% saveas(figure(1), "/Users/leonardo/Desktop/Politecnico/5 anno/Nonlinear Optics/Simulation/Figure1_Ar.png")
% saveas(figure(2), "/Users/leonardo/Desktop/Politecnico/5 anno/Nonlinear Optics/Simulation/Figure2_Ar.png")
% saveas(figure(3), "/Users/leonardo/Desktop/Politecnico/5 anno/Nonlinear Optics/Simulation/Figure3_Ar.png")
% saveas(figure(4), "/Users/leonardo/Desktop/Politecnico/5 anno/Nonlinear Optics/Simulation/Figure4_Ar.png")
% saveas(figure(5), "/Users/leonardo/Desktop/Politecnico/5 anno/Nonlinear Optics/Simulation/Figure5_Ar.png")
% saveas(figure(6), "/Users/leonardo/Desktop/Politecnico/5 anno/Nonlinear Optics/Simulation/Figure6_Ar.png")






%% FUNCTIONS
% Sellmeier equation:
% c = coefficients [dimensionless/microm^2/microm^-2], l = wavelength [microm]
function n = n0(c, l)
    n = zeros(size(l));
    for i=2:2:length(c)
        n = n + c(i)./(1-c(i+1)./l.^2);
    end
    n = n + c(1); % minimize floating point algebra errors by adding the
    % biggest value (in abs) as last one
    n = sqrt(n); % sqrt since Sellmeier eq. gives n^2(lambda)
end

% Second derivative of the refractive index from Sellmeier equation:
% d^2n/dlambda^2 = E(lambda)/n(lambda) + D^2(lambda)/n^3(lambda) [microm^-2]
% D(lambda) = -sum_i{Bi*Ci/(lambda^3*(1-Ci/lambda^2)^2)}
% E(lambda) = -sum_i{4*Bi*Ci^2/(lambda^6*(1-Ci/lambda^2)^3 + 3*Bi*Ci/(lambda^4*(1-Ci/lambda^2)^2)))}
% functions calculated with WolframAlpha
function d2n = d2n0(c, l)
    d2n = E(c, l)./n0(c, l) + (D(c, l).^2)./n0(c, l).^3;
end

function d = D(c, l)
    d = zeros(size(l));
    for i=2:2:length(c)
        d = d + c(i)*c(i+1)./(l.^3.*(1-c(i+1)./l.^2).^2);
    end
end

function e = E(c, l)
    e = zeros(size(l));
    for i=2:2:length(c)
        e = e + 4*c(i)*c(i+1)^2./(l.^6.*(1-c(i+1)./l.^2).^3) + 3*c(i)*c(i+1)./(l.^4.*(1-c(i+1)./l.^2).^2);
    end
end

% Second derivative of the wavevector:
% d^2k0/dw^2|w0 = lambda0^3/(2*pi*c^2)*d^2n/dlambda^2|lambda0 [microm^-1*fs^2]
% where w0 = 2*pi*c/lambda0.
% with [d^2n/dlambda^2] = microm^-2, [lambda0] = microm, [c] = microm/fs
function k = d2k0(c, l)
    v = 299792458e-9; % ~3e8 speed of light [m/s] => *10^6 to [microm/s] => *10^-15 to [microm/fs]
    k = l.^3/(2*pi*v^2).*d2n0(c, l);
end

% Dispersion length:
% [l] = microm, [dt] = fs
function ld = Ld(c, l, dt)
    sigmat = dt/sqrt(8*log(2));
    ld = sigmat^2/abs(d2k0(c, l)); % [microm]
    ld = ld*1e-6; % [m]
end

% Nonlinear length:
% [l] = microm, [E0] = V/m
function ln = Ln(n_0, n_2, l, E0)
    Z0 = 377; % [Ohm] vacuum impedance
    Gbar = 4*pi*n_0*n_2/(Z0*l*1e-6); % [m/V^2]
    ln = 1/(Gbar*E0*conj(E0)); % [m]
end

% Fitting model for Spectral Phase:
% in order to understand how to compensate the spectral phase introduced by
% the dispersion and nonlinear effects, we express it in terms of powers of
% the frequency as phi(w) = a0 + a1*(w-w0) + a2*(w-w0)^2 + a3*(w-w0)^3 + ...
% w & w0 = column frequency vectors, N = order of polynomial
% coeff = column coefficients vector (w^0, w^1, w^2, w^3, ...)
function coeff = PhaseFit(w, w0, N, spectral_phase, spectrum)
    indeces_dw = PulseFind(spectrum);
    w_new = w(indeces_dw(1):indeces_dw(2));
    M = zeros(length(w_new), N);
    M(:,1) = 1;
    for i=2:N
        M(:,i) = (w_new-w0).^(i-1);
    end
    % Ax=b overdetermined system to find coefficients (A=M, x=coeff, b=spectral_phase)
    coeff = M\spectral_phase(indeces_dw(1):indeces_dw(2));
    
end

% Finding the indeces of the vector for which the pulse is <1/1000 of its
% maximum value:
function ind = PulseFind(Ipulse)
    m = max(abs(Ipulse));
    ind = zeros(2,1);
    ind(1) = find((abs(Ipulse) > (1e-3*m)) == 1, 1, 'first');
    ind(2) = find((abs(Ipulse) > (1e-3*m)) == 1, 1, 'last');
end

% Building the Fitted Spectral Phase:
% fit(w) = a0 + a1*(w-w0) + a2*(w-w0)^2 + a3*(w-w0)^3 + ...
function fit = BuildPhaseFit(w, w0, coeff)
    fit = ones(size(w)).*coeff(1);
    for i=2:length(coeff)
        fit = fit + coeff(i).*(w-w0).^(i-1);
    end
    fit = fit';
end

% Dispersion compensation:
% compensate spectral dispersion only in the region where there is actually
% the pulse
function Ecompensated = DispersionCompensation(S_w, phase)
    Ecompensated = fftshift(ifft(S_w.*exp(1i*phase)));
end