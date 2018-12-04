%22.05 F18
%Jacob Miske
%P-Set 3
clc; clear all
%% Problem 1
% Using only the 36.68 eV resonance and the energy range
% [25.0 - 50.0] eV

%compute numerically the U-238 group capture cross section

%Temperature of U-238 and Atomic Mass
A = 238; T=0; %amus, degrees K

%SLBW, equation constants
Eneutron=1000000; %starting energy, eV (1e6 = MeV)
E_1 = 25.0; E_2 = 50.0; E_0 = 36.68; %eV
r_0=(2603911/E_0)*((A+1)/A)^2;
gammaN = 0.03355; gammaY = 0.02300; %eV
gammaTotal = gammaN+gammaY; %eV
%Psi Function assumption
Chi = 2*(Eneutron-E_0)/gammaTotal; Chi1 = 2*(E_1-E_0)/gammaTotal; Chi2 = 2*(E_2-E_0)/gammaTotal;
fChi = @(E) 2*(E-E_0)/gammaTotal;
Psi = (1/(1+(Chi).^2));
fPsi = @(E) (1/(1+(fChi).^2));
%Cross section functions as defined in SLBW
funcsigmaGamma = @(E) sqrt(E_0/E)*r_0*(gammaN/gammaTotal)*(gammaY/gammaTotal)*(1/(1+(2*(E-E_0)/gammaTotal).^2));
sigmaGamma = sqrt(E_0/Eneutron)*r_0*(gammaN/gammaTotal)*(gammaY/gammaTotal)*Psi;

%From analytical, constant 1/E assumption solution of the group capture
%cross section function
sigmaGammaA =(1/E_0)*r_0*((gammaN*gammaY)/(gammaTotal^2))*(gammaTotal/2)*((atan(Chi2)) - (atan(Chi1)))/(log(E_2/E_1));


% Problem 1, Part 3
funcsigmaGamma1overE = @(E) sqrt(E_0./E).*(1./E).*r_0.*(gammaN./gammaTotal).*(gammaY./gammaTotal).*(1./(1+(2.*(E-E_0)./gammaTotal).^2));
func1overE = @(E) 1./E;
sigma_groupCapture = integral(funcsigmaGamma1overE, 25.0, 50.0)./integral(func1overE, 25.0, 50.0);


%% Problem 2
%Including all three resonances contribution to cross sections at all
%energies

%P4
T=0; A=238;
r_0=(2603911/E_0)*((A+1)/A)^2;
%Constants
E_0res1 = 6.67; E_0res2 = 20.87; E_0res3 = 36.68; %eV
E_1res1 = 1.0; E_2res1 = 6.0; %eV
E_1res2 = 6.0; E_2res2 = 10.0; %eV
E_1res3 = 10.0; E_2res3 = 25.0; %eV
E_1res4 =25.0;E_2res4 = 50.0; %eV
gammaNres1 = 0.00148; gammaYres1 = 0.02300; %eV
gammaNres2 = 0.01009; gammaYres2= 0.02286; %eV
gammaNres3= 0.03355; gammaYres3 = 0.02300; %eV
gammaTotalres1 = gammaNres1+gammaYres1; %eV
gammaTotalres2 = gammaNres2+gammaYres2; %eV
gammaTotalres3 = gammaNres3+gammaYres3; %eV
%Three functions, one for each resonance's capture cross section
r_0res1=(2603911/E_0res1)*((A+1)/A)^2;
r_0res2=(2603911/E_0res2)*((A+1)/A)^2;
r_0res3=(2603911/E_0res3)*((A+1)/A)^2;
%define group capture functions for each resonance
funcsigmaGammaRes1times1overE = @(E) sqrt(E_0res1./E).*(1./E).*r_0res1.*(gammaNres1./gammaTotalres1).*(gammaYres1./gammaTotalres1).*(1./(1+(2.*(E-E_0res1)./gammaTotalres1).^2));
funcsigmaGammaRes2times1overE = @(E) sqrt(E_0res2./E).*(1./E).*r_0res2.*(gammaNres2./gammaTotalres2).*(gammaYres2./gammaTotalres2).*(1./(1+(2.*(E-E_0res2)./gammaTotalres2).^2));
funcsigmaGammaRes3times1overE = @(E) sqrt(E_0res3./E).*(1./E).*r_0res3.*(gammaNres3./gammaTotalres3).*(gammaYres3./gammaTotalres3).*(1./(1+(2.*(E-E_0res3)./gammaTotalres3).^2));

func1overE = @(E) 1./E;
sigma_groupCapture = integral(funcsigmaGamma1overE, 25.0, 50.0)./integral(func1overE, 25.0, 50.0);
%Take function from lowest E of 0.1 eV to 50.0 eV
funcsigmaGammaCaptureXSRes1times1overE = @(E) (sqrt(E_0res1./E).*(1./E).*r_0res1.*(gammaNres1./gammaTotalres1).*(gammaYres1./gammaTotalres1).*(1./(1+(2.*(E-E_0res1)./gammaTotalres1).^2))) ./ (1./E);
funcsigmaGammaCaptureXSRes2times1overE = @(E) (sqrt(E_0res2./E).*(1./E).*r_0res2.*(gammaNres2./gammaTotalres2).*(gammaYres2./gammaTotalres2).*(1./(1+(2.*(E-E_0res2)./gammaTotalres2).^2))) ./ (1./E);
funcsigmaGammaCaptureXSRes3times1overE = @(E) (sqrt(E_0res3./E).*(1./E).*r_0res3.*(gammaNres3./gammaTotalres3).*(gammaYres3./gammaTotalres3).*(1./(1+(2.*(E-E_0res3)./gammaTotalres3).^2))) ./ (1./E);
%Take integral from lowest E of 0.1 eV to 50.0 eV
funcsigmaGammaCaptureXSintegralRes1times1overE = @(E) (sqrt(E_0res1./E).*(1./E).*r_0res1.*(gammaNres1./gammaTotalres1).*(gammaYres1./gammaTotalres1).*(1./(1+(2.*(E-E_0res1)./gammaTotalres1).^2))) ./ (1./E);
funcsigmaGammaCaptureXSintegralRes2times1overE = @(E) (sqrt(E_0res2./E).*(1./E).*r_0res2.*(gammaNres2./gammaTotalres2).*(gammaYres2./gammaTotalres2).*(1./(1+(2.*(E-E_0res2)./gammaTotalres2).^2))) ./ (1./E);
funcsigmaGammaCaptureXSintegralRes3times1overE = @(E) (sqrt(E_0res3./E).*(1./E).*r_0res3.*(gammaNres3./gammaTotalres3).*(gammaYres3./gammaTotalres3).*(1./(1+(2.*(E-E_0res3)./gammaTotalres3).^2))) ./ (1./E);


dE = 0.01; numOfSplits = [0.1,dE,50.0]%eV
for 

%plot numerator functions before determining max value of barns
figure(1)
fplot(funcsigmaGammaRes1times1overE, 'r'); hold on
fplot(funcsigmaGammaRes2times1overE, 'b'); hold on
fplot(funcsigmaGammaRes3times1overE, 'g'); grid on
title('Group Capture Cross Section Numerators')
legend('Resonance 1', 'Resonance 2', 'Resonance 3')
xlabel('E (eV)'); ylabel('XS (barns)')
saveas(gcf,'Group Capture Cross Section Numerators.pdf')

%Plot each resonance's impact individually
figure(2)
fplot(funcsigmaGammaCaptureXSRes1times1overE, 'r'); hold on
fplot(funcsigmaGammaCaptureXSRes2times1overE, 'b'); hold on
fplot(funcsigmaGammaCaptureXSRes3times1overE, 'g');
title('Group Capture Cross Section Equation')
legend('Resonance 1', 'Resonance 2', 'Resonance 3')
xlim([0 50]); ylim([0 100]) ;  grid on
xlabel('E (eV)'); ylabel('XS (barns)')
saveas(gcf,'Group Capture Cross Section Equation.pdf')

figure(3)
fplot(funcsigmaGammaCaptureXSintegralRes1times1overE, 'r'); hold on
fplot(funcsigmaGammaCaptureXSintegralRes2times1overE, 'b'); hold on
fplot(funcsigmaGammaCaptureXSintegralRes3times1overE, 'g');
title('Group Capture Cross Section Equation')
legend('Resonance 1', 'Resonance 2', 'Resonance 3')
xlim([0 50]); ylim([0 100]) ;  grid on
xlabel('E (eV)'); ylabel('XS (barns)')
saveas(gcf,'Group Capture Cross Section Equation.pdf')


%Problem 4, Part 2 (Part 5 of pset?)
%compute numerically the resonance integrals for energy groups

%RI_inf^(E1,E2) = int from E1 to E2 of { sigma_gamma(E)  1/E dE

E_1res1 = 1.0; E_2res1 = 6.0; E_1res2 = 6.0; E_2res2 = 10.0; %eV
E_1res3 = 10.0; E_2res3 = 25.0; E_1res4 =25.0;E_2res4 = 50.0; %eV


figure(4)
fplot(funcsigmaGammaCaptureXSRes1times1overE, 'r'); hold on
fplot(funcsigmaGammaCaptureXSRes2times1overE, 'b'); hold on
fplot(funcsigmaGammaCaptureXSRes3times1overE, 'g');
title('Group Capture Cross Section Equation')
legend('Resonance 1', 'Resonance 2', 'Resonance 3')
xlim([0 50]); ylim([0 100]) ;  grid on
xlabel('E (eV)'); ylabel('XS (barns)')
saveas(gcf,'Group Capture Cross Section Equation.pdf')
