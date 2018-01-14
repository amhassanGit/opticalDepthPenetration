%[lam_h2o mua_h2o]=textread('Hale_Querry_water_abs.txt','','headerlines',2);
%[lam_h2o mua_h2o]=textread('wieliczka89.dat','','headerlines',6);
[lam_h2o mua_h2o]=textread('kou93b.dat','','headerlines',6);
lam_h2o = flip(lam_h2o);
mua_h2o = flip(mua_h2o);
%mua_h20 is absorption coefficient (1/cm) of water across wavelength lam_h20

lam=670:5:2400;

% absorption coefficient
water_conc = 0.75; % assume 75% water in brain
mua_water = interp1(lam_h2o, mua_h2o, lam);
% mua is interpolated values of water absorption coefficience at lam wls

% scale mua by water conc and convert to 1/mm by dividing by 10
mua_brain = mua_water/10*water_conc;
mua_brain_Allwater=mua_water/10;
% scattering coefficient
a = 1.1; a_min=1.09; a_max=4.08;% units are 1/mm
b=1.37; b_min=0.334; b_max=3.089; %no units
g = 0.9; % anisotropy
mus_prime = a*(lam/500).^(-b);
mus_prime_min=a_min*(lam/500).^(-b_min);
mus_prime_max=a_max*(lam/500).^(-b_max);
mus = mus_prime/(1-g);
mus_min = mus_prime_min/(1-g);
mus_max = mus_prime_max/(1-g);

%T = exp(-(mua_brain+mus)*z); %Transmission (surviving fraction of light at z)

% Figure for paper
% Fraction of photons arriving at focal volume, normalized to 1300 nm
z= 1; % 1mm
T = exp(-(mua_brain+mus)*z); %Transmission (surviving fraction of light at z for scattering and absorption)
T_min = exp(-(mua_brain+mus_min)*z); %Transmission min (surviving fraction of light at z for scattering and absorption)
T_max = exp(-(mua_brain+mus_max)*z); %Transmission max (surviving fraction of light at z for scattering and absorption)

%Photons absorbed
lam_2=670:0.5:2400;
mua_water_2 = interp1(lam_h2o, mua_h2o, lam_2);
mua_brain_2 = mua_water_2/10*water_conc;
T_abs = exp(-(mua_brain_2)*z); %Transmission (surviving fraction of light at z for absorption)
level=0.5;
aboveLine = (T_abs>=level);
below=T_abs; above=T_abs;
below(aboveLine) = NaN;
above(~aboveLine) = NaN;

close all;
yyaxis left
plot(lam, T/T(127),'b','LineWidth',2); %normalized to 1300 nm (index 127 in lam)
%hold on; plot(lam, T_min/T(127),'b--');
%hold on; plot(lam, T_max/T(127),'b--');
xlabel('Wavelength (nm)') % x-axis label
ylabel('Normalized photon fraction at z=1 mm') % y-axis label
set(gca, 'YColor', [0 0 1]);
hold on; yyaxis right;
plot(lam_2, (1-below)*100,'r',lam_2, (1-above)*100,'r:');
ylabel('Percent photons absorbed');
set(gca, 'YColor', [1 0 0]);
hold off
xlim([700,2400])
ylim([0 102])
set(gca,'fontsize',18)
%legend('','')
%legend('Location','northwest')
legend('boxoff')
ax = gca;
ax.XAxis.MinorTick = 'on';
ax.XAxis.MinorTickValues = ax.XAxis.Limits(1):100:ax.XAxis.Limits(2);

%% Scattering and Absorption Attenuation lengths
%mimic of xu plot
eff=1000./((mus)+(mua_brain));
%close all;
% Plot of scattering length (um) in brain
figure, plot(lam,1./mus*1000,'r');
% Plot of absorption length (um) in brain
hold on; plot(lam,1./mua_brain*1000,'b'); 
% Plot effective attenuation length
hold on; plot(lam, eff, 'g');
% Plot of absorption length (um) in brain if assume 100% water
%hold on; plot(lam,1./mua_brain_Allwater*1000); 
hold off
xlabel('Wavelength (nm)') % x-axis label
ylabel('Attenuation length in brain [um]') % y-axis label
ylim([0,1600])
xlim([700,2000])
%% Transmisssion (surving fraction of light at z)
close all;
zdepth = 0+(0:size(lam,2)-1)*0.002891; % 1 mm depth
T=zeros(347);
for z=1:347
    T(z,:) = exp(-(mua_brain(z)+mus(z))*zdepth);%Transmission (surviving fraction of light at z)
end

figure, surf(zdepth(100:347), lam(100:347), T(100:347,100:347));
ylabel('Wavelength (nm)') % x-axis label
xlabel('z depth (mm)') % y-axis label
zlabel('Transmission at z');
colorbar;

%% Transmisssion (surving fraction of light at z) at single z depth
z= 1; % 1mm
T = exp(-(mua_brain+mus)*z); %Transmission (surviving fraction of light at z)
plot(lam, T);
xlabel('Wavelength (nm)') % x-axis label
ylabel('Transmission at 1 mm') % y-axis label

%% Scattering and Absorption Attenuation coefficients
%Plot of scattering coefficient (1/mm) in brain
close all;
figure, plot(lam,mus); %Plot scattering coefficient
hold on; plot(lam,mua_brain); %Plot absorption coefficient
hold on; plot(lam,1./((1./mus)+(1./mua_brain)), 'g'); %Plot effective attenuation coefficient
xlabel('Wavelength (nm)') % x-axis label
ylabel('Attenuation coefficents in brain (1/mm)') % y-axis label

%% Scattering
%Plot of scattering coefficient (1/mm) in brain
close all;
figure, plot(lam,mus);
xlabel('Wavelength (nm)') % x-axis label
ylabel('Scattering coefficent in brain (1/mm)') % y-axis label

% Plot of scattering length (um) in brain
figure, plot(lam,1./mus*1000);
xlabel('Wavelength (nm)') % x-axis label
ylabel('Scattering length in brain (um)') % y-axis label
ylim([0,2000])
xlim([600,2000])

%% Absorption 
%Plot of absorption coefficient (1/mm) in brain
close all;
figure, plot(lam,mua_brain);
xlabel('Wavelength (nm)') % x-axis label
ylabel('Absorption coefficent in brain (1/mm)') % y-axis label

% Plot of absorption length (um) in brain
figure, plot(lam,1./mua_brain*1000);
xlabel('Wavelength (nm)') % x-axis label
ylabel('Absorption length for water in brain (um)') % y-axis label
ylim([0,2000])
xlim([600,2000])