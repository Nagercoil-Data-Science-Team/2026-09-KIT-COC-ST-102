clc; 
clear; 
close all;

% Wider default figure size so titles/labels are never clipped on export
set(0,'DefaultFigurePosition',[100 100 1400 900]);

%% ============================================================
%% ============================================================
%% ============================================================
% STEP 1 — A/B-SITE DOPING MATRIX
%% ============================================================

x_values = 0:0.05:0.70;
y_values = 0:0.05:0.70;
[X,Y] = meshgrid(x_values,y_values);


Vo_component = 0.10 + 0.30*X + 0.25*Y;                 % oxygen-vacancy proxy
OEC_component = 0.08 + 0.35*X + 0.30*Y;                % oxygen-exchange-capacity proxy
EnergyDensity_component = 50 + 100*X + 90*Y;           % thermochemical storage / energy-density proxy
Stability_component = 90 - 10*X - 8*Y;                 % structural/redox stability proxy

w_vacancy = 0.30;
w_oec = 0.30;
w_energy = 0.30;
w_stability = 0.10;

Vo_norm = (Vo_component-min(Vo_component(:)))/(max(Vo_component(:))-min(Vo_component(:)));
OEC_norm = (OEC_component-min(OEC_component(:)))/(max(OEC_component(:))-min(OEC_component(:)));
Energy_norm = (EnergyDensity_component-min(EnergyDensity_component(:)))/(max(EnergyDensity_component(:))-min(EnergyDensity_component(:)));
Stability_norm = (Stability_component-min(Stability_component(:)))/(max(Stability_component(:))-min(Stability_component(:)));

optimization_score = w_vacancy*Vo_norm + w_oec*OEC_norm + w_energy*Energy_norm + w_stability*Stability_norm;

% Kept for backward-compatible reporting/plotting under the original names
oxygen_vacancy_score = Vo_component;
oxygen_exchange_score = OEC_component;
storage_score = EnergyDensity_component;
stability_score = Stability_component;

%% Find optimum

[best_score,best_index] = max(optimization_score(:));
best_x = X(best_index);
best_y = Y(best_index);

%% Formulation

best_formula = sprintf('La_{%.2f}Ca_{%.2f}Co_{%.2f}Fe_{%.2f}O_3',1-best_x,best_x,1-best_y,best_y);

fprintf('\n====================================================\n');
fprintf('              OPTIMIZATION RESULT\n');
fprintf('====================================================\n');
fprintf('Weights: w_vacancy=%.2f, w_OEC=%.2f, w_energy=%.2f, w_stability=%.2f (sum=%.2f)\n',w_vacancy,w_oec,w_energy,w_stability,w_vacancy+w_oec+w_energy+w_stability);
fprintf('NOTE: All four objectives are min-max normalized before weighting.\n');
fprintf('Grid tested: x = %.2f:%.2f:%.2f, y = %.2f:%.2f:%.2f\n',min(x_values),x_values(2)-x_values(1),max(x_values),min(y_values),y_values(2)-y_values(1),max(y_values));
fprintf('Best x value            = %.2f\n',best_x);
fprintf('Best y value            = %.2f\n',best_y);
fprintf('Best optimization score = %.6f (0-1 normalized scale under this model/weighting; not a percentage performance)\n',best_score);
fprintf('\nBest formulation (optimal only under this assumed model, not yet experimentally validated):\n');
fprintf('La_{%.2f}Ca_{%.2f}Co_{%.2f}Fe_{%.2f}O3\n',1-best_x,best_x,1-best_y,best_y);
if best_x == max(x_values) && best_y == max(y_values)
    fprintf('\nNOTE: Optimum still lies at the upper boundary of the tested x/y range (%.2f).\n',max(x_values));
    fprintf('Extend the search grid further to confirm whether this is a true optimum\n');
    fprintf('or an artifact of the tested range.\n');
end
fprintf('\nOPTIMIZATION CRITERION (Score A): weighted, min-max-normalized combination of\n');
fprintf('oxygen-vacancy proxy, OEC proxy, energy-density proxy, and stability proxy\n');
fprintf('(weights %.2f/%.2f/%.2f/%.2f). This is a DIFFERENT scoring criterion from the\n',w_vacancy,w_oec,w_energy,w_stability);
fprintf('Structure-Vacancy-Storage combined score computed later (Step 27) — the two are\n');
fprintf('not expected to agree, since they weight different and only partially overlapping\n');
fprintf('quantities. See the OPTIMIZATION CRITERIA COMPARISON section at the end of this\n');
fprintf('run for a side-by-side explanation.\n');


%% ============================================================
% STEP 3 — SIMULATED XRD ANALYSIS
%% ============================================================

fprintf('\n====================================================\n');
fprintf('                  XRD ANALYSIS\n');
fprintf('====================================================\n');

two_theta = 20:0.02:80;
peak_positions = [32.50 46.70 58.20 68.50];
peak_intensities = [100 65 45 30];
peak_width = 0.24;
xrd_intensity = zeros(size(two_theta));

for i = 1:length(peak_positions)
    xrd_intensity = xrd_intensity + peak_intensities(i) .* exp(-4*log(2) .* ((two_theta-peak_positions(i))./peak_width).^2);
end

xrd_intensity = xrd_intensity + 3;

%% Scherrer equation

K = 0.9;
lambda = 0.15406;        % nm, Cu K-alpha
theta = deg2rad(peak_positions(1)/2);
beta_rad = deg2rad(peak_width);
crystallite_size = K*lambda/(beta_rad*cos(theta));

fprintf('Main peak position     = %.2f degree\n',peak_positions(1));
fprintf('Peak width              = %.4f degree\n',peak_width);
fprintf('Crystallite size        = %.3f nm\n',crystallite_size);

%% XRD plot

fig = figure('Name','XRD - Best Formulation');
plot(two_theta,xrd_intensity,'LineWidth',2,'Color',[0.85 0.10 0.10]);
xlabel('2\theta (degree)','FontName','Times New Roman','FontSize',20,'FontWeight','bold');
ylabel('Intensity (a.u.)','FontName','Times New Roman','FontSize',20,'FontWeight','bold');
title('XRD - Best Formulation','FontName','Times New Roman','FontSize',20,'FontWeight','bold');
ax = gca; ax.FontName = 'Times New Roman'; ax.FontSize = 20; ax.FontWeight = 'bold';
drawnow;
exportgraphics(fig,'XRD_best_formulation.png','Resolution',1000);


%% ============================================================
% STEP 4 — SIMULATED SEM ANALYSIS
%% ============================================================

fprintf('\n====================================================\n');
fprintf('                  SEM ANALYSIS\n');
fprintf('====================================================\n');

rng(10);
number_particles = 250;
particle_size = 105 + 18*randn(number_particles,1);
particle_size = max(particle_size,50);
mean_particle_size = mean(particle_size);
median_particle_size = median(particle_size);
std_particle_size = std(particle_size);
minimum_particle_size = min(particle_size);
maximum_particle_size = max(particle_size);

fprintf('Number of particles    = %d\n',number_particles);
fprintf('Mean particle size     = %.3f nm\n',mean_particle_size);
fprintf('Median particle size   = %.3f nm\n',median_particle_size);
fprintf('Std. deviation         = %.3f nm\n',std_particle_size);
fprintf('Minimum particle size  = %.3f nm\n',minimum_particle_size);
fprintf('Maximum particle size  = %.3f nm\n',maximum_particle_size);

%% SEM particle-size plot

fig = figure('Name','SEM Particle Size');
histogram(particle_size,15,'FaceColor',[0.20 0.60 0.20],'EdgeColor',[0 0 0]);
xlabel('Particle size (nm)','FontName','Times New Roman','FontSize',20,'FontWeight','bold');
ylabel('Frequency','FontName','Times New Roman','FontSize',20,'FontWeight','bold');
title('SEM Particle-Size Distribution','FontName','Times New Roman','FontSize',20,'FontWeight','bold');
ax = gca; ax.FontName = 'Times New Roman'; ax.FontSize = 20; ax.FontWeight = 'bold';
drawnow;
exportgraphics(fig,'SEM_particle_size.png','Resolution',1000);

%% Simulated SEM morphology

image_size = 500;
sem_image = zeros(image_size);

for i = 1:120
    cx = randi([20 image_size-20]);
    cy = randi([20 image_size-20]);
    radius = randi([5 18]);
    [xx,yy] = meshgrid(1:image_size,1:image_size);
    particle = exp(-((xx-cx).^2 + (yy-cy).^2)/(2*radius^2));
    sem_image = sem_image + particle;
end

sem_image = sem_image ./ max(sem_image(:));

fig = figure('Name','SEM Morphology');
imagesc(sem_image); axis image off; colormap(gca,'parula');
title('SEM Morphology - Best Formulation','FontName','Times New Roman','FontSize',20,'FontWeight','bold');
ax = gca; ax.FontName = 'Times New Roman'; ax.FontSize = 20; ax.FontWeight = 'bold';
drawnow;
exportgraphics(fig,'SEM_morphology.png','Resolution',1000);


%% ============================================================
% STEP 5 — SIMULATED BET ANALYSIS
%% ============================================================

fprintf('\n====================================================\n');
fprintf('                  BET ANALYSIS\n');
fprintf('====================================================\n');

relative_pressure = (0.01:0.01:0.99)';
adsorption_volume = 5 + 18*relative_pressure ./ (1-relative_pressure);
desorption_volume = 5 + 16*relative_pressure ./ (1-relative_pressure);
BET_surface_area = 32.0;
total_pore_volume = 0.190;
average_pore_diameter = 4*total_pore_volume*1000/BET_surface_area;

fprintf('BET surface area       = %.3f m^2/g\n',BET_surface_area);
fprintf('Total pore volume      = %.5f cm^3/g\n',total_pore_volume);
fprintf('Average pore diameter  = %.3f nm\n',average_pore_diameter);

%% BET plot

fig = figure('Name','BET Isotherm');
plot(relative_pressure,adsorption_volume,'LineWidth',2,'Color',[0.00 0.45 0.74]);
hold on;
plot(relative_pressure,desorption_volume,'LineWidth',2,'Color',[0.93 0.49 0.19]);
hold off;
xlabel('Relative pressure (P/P_0)','FontName','Times New Roman','FontSize',20,'FontWeight','bold');
ylabel('Adsorbed volume (a.u.)','FontName','Times New Roman','FontSize',20,'FontWeight','bold');
title('BET Adsorption-Desorption Isotherm','FontName','Times New Roman','FontSize',20,'FontWeight','bold');
legend({'Adsorption','Desorption'},'FontName','Times New Roman','FontSize',16,'FontWeight','bold');
ax = gca; ax.FontName = 'Times New Roman'; ax.FontSize = 20; ax.FontWeight = 'bold';
drawnow;
exportgraphics(fig,'BET_isotherm.png','Resolution',1000);


%% ============================================================
% STEP 6 — SIMULATED/REFERENCE DELTA0 (NOT A REAL IODOMETRIC RESULT)

%% ============================================================

fprintf('\n====================================================\n');
fprintf('          OXYGEN NON-STOICHIOMETRY\n');
fprintf('====================================================\n');

M_La = 138.90547; M_Ca = 40.078; M_Co = 58.933; M_Fe = 55.845; M_O = 15.999;
M_formula = (1-best_x)*M_La + best_x*M_Ca + (1-best_y)*M_Co + best_y*M_Fe + 3*M_O;

sample_mass = 0.05000;       % g
titrant_volume = 12.50;      % mL
titrant_molarity = 0.0100;   % mol/L

n_sample = sample_mass/M_formula;
n_titrant = (titrant_volume/1000)*titrant_molarity;
n_I2 = n_titrant/2;
delta_0 = n_I2/n_sample;
delta_0 = max(0,min(delta_0,0.50));

fprintf('\n--- SIMULATED/REFERENCE DELTA0 (placeholder, not a real iodometric titration) ---\n');
fprintf('Sample mass             = %.5f g\n',sample_mass);
fprintf('Titrant volume          = %.3f mL\n',titrant_volume);
fprintf('Titrant concentration   = %.4f mol/L\n',titrant_molarity);
fprintf('Formula molar mass      = %.4f g/mol\n',M_formula);
fprintf('Simulated/Reference delta0 = %.5f\n',delta_0);


%% ============================================================
% STEP 7 — CONSISTENT TGA DELTA MODEL (REALISTIC RELEASE SHAPE)

%% ============================================================

temperature = (300:5:1000)';
Delta_delta_best = 0.10 + 0.10*best_x + 0.08*best_y;
Delta_delta_best = min(Delta_delta_best,0.20);

% Three overlapping release contributions (illustrative TPD shape):
%   low-T release  ~400-600 C, main release ~600-850 C, high-T tail ~800-950 C
release_peak_centers = [500 720 880];
release_peak_widths = [55 60 55];
release_peak_weights = [0.20 0.50 0.30];   % must sum to 1

release_rate_shape = zeros(size(temperature));
for k = 1:numel(release_peak_centers)
    release_rate_shape = release_rate_shape + release_peak_weights(k)*exp(-4*log(2)*((temperature-release_peak_centers(k))/release_peak_widths(k)).^2);
end

release_rate = (release_rate_shape/trapz(temperature,release_rate_shape))*Delta_delta_best;
Delta_delta_T = cumtrapz(temperature,release_rate);
delta_absolute = delta_0 + Delta_delta_T;
vacancy_fraction = delta_absolute/3;
vacancy_percentage = vacancy_fraction*100;
oxygen_mass_fraction = Delta_delta_T*M_O/M_formula;
mass_TGA = sample_mass .* (1-oxygen_mass_fraction);

fprintf('\n--- TGA OXYGEN NON-STOICHIOMETRY ---\n');
fprintf('Initial delta0          = %.5f\n',delta_0);
fprintf('Maximum Delta-delta     = %.5f\n',Delta_delta_best);
fprintf('Final absolute delta    = %.5f\n',delta_absolute(end));
fprintf('Final vacancy fraction  = %.5f (= delta/3, three O sites per ABO3 formula unit)\n',vacancy_fraction(end));
fprintf('Final oxygen-site deficiency = %.3f %% of the three nominal oxygen sites\n',vacancy_percentage(end));

%% TGA / Delta plot

fig = figure('Name','Oxygen Non-Stoichiometry');
plot(temperature,delta_absolute,'LineWidth',2,'Color',[0.49 0.18 0.56]);
xlabel('Temperature (°C)','FontName','Times New Roman','FontSize',20,'FontWeight','bold');
ylabel('\delta','FontName','Times New Roman','FontSize',20,'FontWeight','bold');
title('Oxygen Non-Stoichiometry of Best Formulation','FontName','Times New Roman','FontSize',20,'FontWeight','bold');
ax = gca; ax.FontName = 'Times New Roman'; ax.FontSize = 20; ax.FontWeight = 'bold';
drawnow;
exportgraphics(fig,'Oxygen_NonStoichiometry.png','Resolution',1000);


%% ============================================================
% STEP 8 — A/B-SITE DOPING VS OXYGEN VACANCIES
% SAME DELTA MODEL AS THE TGA MODEL
%% ============================================================

results_AB = [];

for i = 1:length(x_values)
    for j = 1:length(y_values)
        x_current = x_values(i);
        y_current = y_values(j);
        Delta_delta_current = 0.10 + 0.10*x_current + 0.08*y_current;
        Delta_delta_current = min(Delta_delta_current,0.20);
        delta_current = delta_0 + Delta_delta_current;
        vacancy_current = delta_current/3;
        vacancy_percent_current = vacancy_current*100;
        results_AB = [results_AB; x_current, y_current, Delta_delta_current, delta_current, vacancy_current, vacancy_percent_current];
    end
end

[maximum_vacancy,best_vacancy_index] = max(results_AB(:,6));
vacancy_x = results_AB(best_vacancy_index,1);
vacancy_y = results_AB(best_vacancy_index,2);
vacancy_delta = results_AB(best_vacancy_index,4);

fprintf('\n====================================================\n');
fprintf('      EFFECT OF A/B-SITE DOPING ON VACANCIES\n');
fprintf('====================================================\n');
fprintf('Highest oxygen-site deficiency composition:\n');
fprintf('x = %.2f\n',vacancy_x);
fprintf('y = %.2f\n',vacancy_y);
fprintf('delta = %.5f\n',vacancy_delta);
fprintf('Vo/3  = %.5f  (i.e. %.3f %% of the three nominal oxygen sites)\n',maximum_vacancy/100,maximum_vacancy);
fprintf('Suggested reporting: "The oxygen-site deficiency corresponds to approximately %.2f%% of the three nominal oxygen sites (delta = %.4f)."\n',maximum_vacancy,vacancy_delta);


%% ============================================================
% STEP 9 — CONSISTENCY CHECK
%% ============================================================

best_row = results_AB(results_AB(:,1)==best_x & results_AB(:,2)==best_y,:);
delta_from_AB = best_row(4);
vacancy_from_AB = best_row(6);
delta_from_TGA = delta_absolute(end);
vacancy_from_TGA = vacancy_percentage(end);

fprintf('\n====================================================\n');
fprintf('                 CONSISTENCY CHECK\n');
fprintf('====================================================\n');
fprintf('Best formulation x       = %.2f\n',best_x);
fprintf('Best formulation y       = %.2f\n',best_y);
fprintf('\nTGA absolute delta       = %.5f\n',delta_from_TGA);
fprintf('A/B absolute delta       = %.5f\n',delta_from_AB);
fprintf('\nTGA vacancy               = %.3f %%\n',vacancy_from_TGA);
fprintf('A/B vacancy               = %.3f %%\n',vacancy_from_AB);

if abs(delta_from_TGA-delta_from_AB) < 1e-10
    fprintf('\nSTATUS: CONSISTENT\n');
else
    fprintf('\nSTATUS: INCONSISTENT\n');
end


%% ============================================================
% STEP 10 — CO/FE OXIDATION STATES (SIMULATED XPS ANALYSIS)
%% ============================================================

fprintf('\n====================================================\n');
fprintf('    Co / Fe OXIDATION STATE ANALYSIS (SIMULATED XPS)\n');
fprintf('====================================================\n');

Co2_fraction = 0.30 + 0.15*best_y;
Co3_fraction = 1-Co2_fraction;
Fe2_fraction = 0.20 + 0.10*best_y;
Fe3_fraction = 1-Fe2_fraction;
Co_average_valence = 2*Co2_fraction + 3*Co3_fraction;
Fe_average_valence = 2*Fe2_fraction + 3*Fe3_fraction;

fprintf('\nCo oxidation states:\n');
fprintf('Co2+ fraction = %.2f %%\n',Co2_fraction*100);
fprintf('Co3+ fraction = %.2f %%\n',Co3_fraction*100);
fprintf('Average Co valence = %.3f+\n',Co_average_valence);
fprintf('\nFe oxidation states:\n');
fprintf('Fe2+ fraction = %.2f %%\n',Fe2_fraction*100);
fprintf('Fe3+ fraction = %.2f %%\n',Fe3_fraction*100);
fprintf('Average Fe valence = %.3f+\n',Fe_average_valence);


%% ============================================================
% STEP 11 — O 1s XPS (SIMULATED, NOT PEAK-FITTED)
%% ============================================================

lattice_oxygen_fraction = 0.70 - 0.15*(best_x+best_y);
vacancy_related_fraction = 0.15 + 0.10*(best_x+best_y);
adsorbed_oxygen_fraction = 1-lattice_oxygen_fraction-vacancy_related_fraction;

fprintf('\n--- O 1s XPS ANALYSIS (SIMULATED) ---\n');
fprintf('Lattice oxygen       = %.2f %%\n',lattice_oxygen_fraction*100);
fprintf('Vacancy-related O    = %.2f %%\n',vacancy_related_fraction*100);
fprintf('Adsorbed oxygen      = %.2f %%\n',adsorbed_oxygen_fraction*100);

%% O 1s plot

fig = figure('Name','O 1s XPS');
oxygen_species = categorical({'Lattice O','Vacancy-related O','Adsorbed O'});
oxygen_values = [lattice_oxygen_fraction vacancy_related_fraction adsorbed_oxygen_fraction]*100;
b1 = bar(oxygen_species,oxygen_values,'FaceColor','flat');
b1.CData(1,:) = [0.00 0.45 0.74];
b1.CData(2,:) = [0.85 0.33 0.10];
b1.CData(3,:) = [0.47 0.67 0.19];
text(b1.XEndPoints,b1.YEndPoints,string(round(oxygen_values,2))+"%",'HorizontalAlignment','center','VerticalAlignment','bottom','FontName','Times New Roman','FontSize',14,'FontWeight','bold');
ylabel('Relative contribution (%)','FontName','Times New Roman','FontSize',20,'FontWeight','bold');
title('O 1s XPS Oxygen Species','FontName','Times New Roman','FontSize',20,'FontWeight','bold');
ax = gca; ax.FontName = 'Times New Roman'; ax.FontSize = 20; ax.FontWeight = 'bold';
drawnow;
exportgraphics(fig,'O1s_XPS_Oxygen_Species.png','Resolution',1000);


%% ============================================================
% STEP 12 — Co / Fe OXIDATION PLOT
%% ============================================================

fig = figure('Name','Co Fe Oxidation State');
oxidation_values = [Co2_fraction Co3_fraction; Fe2_fraction Fe3_fraction]*100;
b2 = bar(oxidation_values);
b2(1).FaceColor = [0.30 0.75 0.93];
b2(2).FaceColor = [0.64 0.08 0.18];
for k = 1:numel(b2)
    text(b2(k).XEndPoints,b2(k).YEndPoints,string(round(b2(k).YData,2))+"%",'HorizontalAlignment','center','VerticalAlignment','bottom','FontName','Times New Roman','FontSize',13,'FontWeight','bold');
end
xlabel('Transition-metal species','FontName','Times New Roman','FontSize',20,'FontWeight','bold');
ylabel('Relative fraction (%)','FontName','Times New Roman','FontSize',20,'FontWeight','bold');
title('Co and Fe Oxidation States','FontName','Times New Roman','FontSize',20,'FontWeight','bold');
legend({'2+','3+'},'FontName','Times New Roman','FontSize',16,'FontWeight','bold');
set(gca,'XTickLabel',{'Co','Fe'});
ax = gca; ax.FontName = 'Times New Roman'; ax.FontSize = 20; ax.FontWeight = 'bold';
drawnow;
exportgraphics(fig,'Co_Fe_Oxidation_State.png','Resolution',1000);


%% ============================================================
% STEP 13 — A/B DOPING VACANCY SURFACE
%% ============================================================

vacancy_matrix = reshape(results_AB(:,6),length(y_values),length(x_values));
x_tick_step_ab = max(x_values)/7;                  % ~8 tick marks regardless of grid resolution
x_ticks_sparse_ab = 0:x_tick_step_ab:max(x_values);

fig = figure('Name','A/B Doping Oxygen Vacancy');
line_colors = lines(length(y_values));
hold on;
for j = 1:length(y_values)
    plot(x_values,vacancy_matrix(j,:),'-o','LineWidth',1.8,'MarkerSize',4,'Color',line_colors(j,:),'MarkerFaceColor',line_colors(j,:),'DisplayName',sprintf('y = %.2f',y_values(j)));
end
hold off;
xlabel('Ca A-site doping, x','FontName','Times New Roman','FontSize',20,'FontWeight','bold');
ylabel('Oxygen vacancy (%)','FontName','Times New Roman','FontSize',20,'FontWeight','bold');
title('A/B-Site Doping vs Oxygen Vacancy Formation','FontName','Times New Roman','FontSize',18,'FontWeight','bold');
xticks(x_ticks_sparse_ab);
legend('show','Location','eastoutside','NumColumns',2,'FontName','Times New Roman','FontSize',9,'FontWeight','bold');
ax = gca; ax.FontName = 'Times New Roman'; ax.FontSize = 16; ax.FontWeight = 'bold';
drawnow;
exportgraphics(fig,'AB_Doping_Oxygen_Vacancy.png','Resolution',1000);


%% ============================================================
% STEP 14 — DOPING COMPOSITION VS VACANCY
%% ============================================================

fig = figure('Name','Doping vs Oxygen Vacancy');
plot(1:size(results_AB,1),results_AB(:,6),'-o','LineWidth',2,'Color',[0.93 0.69 0.13],'MarkerFaceColor',[0.85 0.33 0.10]);
xlabel('A/B-site doping composition','FontName','Times New Roman','FontSize',20,'FontWeight','bold');
ylabel('Oxygen vacancy (%)','FontName','Times New Roman','FontSize',20,'FontWeight','bold');
title('Effect of A/B-Site Doping on Oxygen Vacancy Formation','FontName','Times New Roman','FontSize',20,'FontWeight','bold');
ax = gca; ax.FontName = 'Times New Roman'; ax.FontSize = 20; ax.FontWeight = 'bold';
drawnow;
exportgraphics(fig,'Doping_vs_Oxygen_Vacancy.png','Resolution',1000);


%% ============================================================
% STEP 15 — OXYGEN RELEASE (TPD-STYLE) PROFILE
% TPD/TGA release -> O2 uptake -> Redox behavior
% Determine oxygen-release temperature and amount
%
% NOTE: Uses the SAME analytic release_rate curve built in Step 7
% (three overlapping Gaussian contributions), so the reported
% onset/peak/end temperatures are read directly off a curve with a
% genuine interior maximum, not from a numerical derivative of a
% linear ramp.
%% ============================================================

fprintf('\n====================================================\n');
fprintf('       OXYGEN RELEASE (TPD/TGA) ANALYSIS\n');
fprintf('====================================================\n');

[max_release_rate,idx_release] = max(release_rate);
T_release_onset = temperature(find(release_rate > 0.05*max_release_rate,1,'first'));
T_release_peak = temperature(idx_release);
T_release_end = temperature(find(release_rate > 0.05*max_release_rate,1,'last'));
total_O2_released_mol_per_g = (Delta_delta_best/2)/M_formula;
total_O2_released_cm3_g = total_O2_released_mol_per_g*22414;
total_O2_released_wtpercent = oxygen_mass_fraction(end)*100;

fprintf('Onset release temperature   = %.1f degC\n',T_release_onset);
fprintf('Peak release temperature    = %.1f degC\n',T_release_peak);
fprintf('End release temperature     = %.1f degC\n',T_release_end);
fprintf('Total O2 released           = %.5f mol O2/g\n',total_O2_released_mol_per_g);
fprintf('Total O2 released           = %.3f cm^3(STP)/g\n',total_O2_released_cm3_g);
fprintf('Total O2 released           = %.4f wt%%\n',total_O2_released_wtpercent);

fig = figure('Name','Oxygen Release TPD Profile');
plot(temperature,release_rate,'LineWidth',2,'Color',[0.85 0.10 0.10]);
hold on;
plot(T_release_peak,max_release_rate,'o','MarkerSize',10,'MarkerFaceColor',[0.00 0.00 0.00],'MarkerEdgeColor',[0.00 0.00 0.00]);
xline(T_release_onset,'--','Onset','Color',[0.4 0.4 0.4],'FontName','Times New Roman','FontSize',12);
xline(T_release_end,'--','End','Color',[0.4 0.4 0.4],'FontName','Times New Roman','FontSize',12);
hold off;
xlabel('Temperature (°C)','FontName','Times New Roman','FontSize',20,'FontWeight','bold');
ylabel('d\delta/dT (release rate, a.u.)','FontName','Times New Roman','FontSize',20,'FontWeight','bold');
title('TPD Oxygen-Release Profile','FontName','Times New Roman','FontSize',20,'FontWeight','bold');
legend({'Release rate','Peak release'},'FontName','Times New Roman','FontSize',14,'FontWeight','bold');
ax = gca; ax.FontName = 'Times New Roman'; ax.FontSize = 20; ax.FontWeight = 'bold';
drawnow;
exportgraphics(fig,'Oxygen_Release_TPD_Profile.png','Resolution',1000);


%% ============================================================
% STEP 16 — REVERSIBLE OXYGEN UPTAKE (RE-OXIDATION)   [FIXED]
% Evaluate reversible oxygen uptake on cooling/re-oxidation
%
% The re-oxidation (cooling) branch mirrors the same release-rate
% SHAPE built in Step 7 (not a separate linear ramp), scaled by a
% doping-dependent reversibility fraction eta_reversibility < 1 to
% represent an irreversible fraction of the released oxygen.
%
% delta_absolute_reox starts at delta_absolute(end) (= delta_max, the
% same point where the heating curve ends) and decreases as
% reox_uptake_fraction grows toward 1.
%% ============================================================

fprintf('\n====================================================\n');
fprintf('       REVERSIBLE OXYGEN UPTAKE (RE-OXIDATION)\n');
fprintf('====================================================\n');

eta_reversibility = 0.90 - 0.05*best_x - 0.03*best_y;
eta_reversibility = max(0,min(eta_reversibility,1));

delta_max = delta_absolute(end);

temperature_cooling = flip(temperature);
cumulative_release_fraction = Delta_delta_T/Delta_delta_best;
reox_uptake_fraction = 1 - flip(cumulative_release_fraction);   % 0 at hot start of cooling -> 1 at full cool

reversible_delta = eta_reversibility*Delta_delta_best;
irreversible_delta = Delta_delta_best - reversible_delta;

Delta_delta_T_reox = reversible_delta*reox_uptake_fraction;     % O2 re-absorbed so far during cooling
delta_absolute_reox = delta_max - Delta_delta_T_reox;            % continuous with heating curve at delta_max

delta_reoxidized_final = delta_absolute_reox(end);
irreversible_delta_loss = delta_reoxidized_final - delta_0;

fprintf('Initial delta0                  = %.5f\n',delta_0);
fprintf('Maximum delta after reduction   = %.5f\n',delta_max);
fprintf('Total Delta-delta               = %.5f\n',Delta_delta_best);
fprintf('\nReversibility fraction          = %.3f\n',eta_reversibility);
fprintf('Reversible Delta-delta          = %.5f\n',reversible_delta);
fprintf('Irreversible Delta-delta        = %.5f\n',irreversible_delta);
fprintf('\nDelta after re-oxidation        = %.5f\n',delta_reoxidized_final);
fprintf('Irreversible delta loss         = %.5f\n',irreversible_delta_loss);
fprintf('Reversible oxygen uptake (%%)     = %.3f %%\n',eta_reversibility*100);

if abs(irreversible_delta_loss-irreversible_delta) < 1e-10 && ...
   abs((delta_max-reversible_delta)-delta_reoxidized_final) < 1e-10
    fprintf('\nSTATUS: INTERNALLY CONSISTENT\n');
else
    fprintf('\nSTATUS: INCONSISTENT\n');
end

fig = figure('Name','Redox Hysteresis - Release and Uptake');
plot(temperature,delta_absolute,'LineWidth',2,'Color',[0.85 0.10 0.10]);
hold on;
plot(temperature_cooling,delta_absolute_reox,'LineWidth',2,'Color',[0.00 0.45 0.74]);
hold off;
xlabel('Temperature (°C)','FontName','Times New Roman','FontSize',20,'FontWeight','bold');
ylabel('\delta','FontName','Times New Roman','FontSize',20,'FontWeight','bold');
title('Oxygen Release (Heating) and Uptake (Cooling)','FontName','Times New Roman','FontSize',20,'FontWeight','bold');
legend({'Release (heating)','Re-oxidation (cooling)'},'FontName','Times New Roman','FontSize',14,'FontWeight','bold');
ax = gca; ax.FontName = 'Times New Roman'; ax.FontSize = 20; ax.FontWeight = 'bold';
drawnow;
exportgraphics(fig,'Redox_Hysteresis_Release_Uptake.png','Resolution',1000);


%% ============================================================
% STEP 17 — OXYGEN EXCHANGE CAPACITY (OEC)
% Calculate oxygen exchange capacity
%
% NOTE: The OEC unit conversions below (mol O2/g <-> mmol/g <->
% cm3(STP)/g <-> wt%) are mutually consistent and do not need to
% change. What SHOULD change once real data exists is the INPUT
% reversible_delta, which here comes from the simulated
% Delta_delta_best * eta_reversibility (now computed once in Step 16
% and reused here). Replace reversible_delta with the experimentally
% measured reversible oxygen exchange from TGA/redox cycling once
% available.
%% ============================================================

fprintf('\n====================================================\n');
fprintf('           OXYGEN EXCHANGE CAPACITY (OEC)\n');
fprintf('====================================================\n');

n_O2_exchange_mol_g = (reversible_delta/2)/M_formula;
OEC_cm3_g = n_O2_exchange_mol_g*22414;
OEC_mmol_g = n_O2_exchange_mol_g*1000;
OEC_wtpercent = (reversible_delta*M_O/M_formula)*100;
irreversible_wtpercent = (irreversible_delta*M_O/M_formula)*100;

fprintf('Reversible Delta-delta (SIMULATED input) = %.5f\n',reversible_delta);
fprintf('Irreversible Delta-delta                 = %.5f\n',irreversible_delta);
fprintf('Predicted oxygen exchange capacity       = %.5f mol O2/g\n',n_O2_exchange_mol_g);
fprintf('Predicted oxygen exchange capacity       = %.4f mmol O2/g\n',OEC_mmol_g);
fprintf('Predicted oxygen exchange capacity       = %.3f cm^3(STP)/g\n',OEC_cm3_g);
fprintf('Predicted oxygen exchange capacity       = %.4f wt%%\n',OEC_wtpercent);
fprintf('NOTE: replace reversible_delta with an experimental TGA/redox-cycling value for real OEC.\n');

fig = figure('Name','Oxygen Release vs Exchange Capacity');
oec_categories = categorical({'Total Released','Reversible (OEC)','Irreversible Loss'});
oec_categories = reordercats(oec_categories,{'Total Released','Reversible (OEC)','Irreversible Loss'});
oec_values = [total_O2_released_wtpercent OEC_wtpercent irreversible_wtpercent];
b3 = bar(oec_categories,oec_values,'FaceColor','flat');
b3.CData(1,:) = [0.47 0.25 0.80];
b3.CData(2,:) = [0.00 0.60 0.30];
b3.CData(3,:) = [0.80 0.20 0.20];
text(b3.XEndPoints,b3.YEndPoints,string(round(oec_values,3))+"%",'HorizontalAlignment','center','VerticalAlignment','bottom','FontName','Times New Roman','FontSize',14,'FontWeight','bold');
ylabel('Oxygen mass fraction (wt%)','FontName','Times New Roman','FontSize',20,'FontWeight','bold');
title('Total O2 Released vs Oxygen Exchange Capacity','FontName','Times New Roman','FontSize',20,'FontWeight','bold');
ax = gca; ax.FontName = 'Times New Roman'; ax.FontSize = 20; ax.FontWeight = 'bold';
drawnow;
exportgraphics(fig,'Oxygen_Exchange_Capacity.png','Resolution',1000);


%% ============================================================
% STEP 18 — DOPING INFLUENCE ON OXYGEN-RELEASE BEHAVIOR
% Investigate the influence of doping on oxygen-release behavior
%% ============================================================

fprintf('\n====================================================\n');
fprintf('     DOPING INFLUENCE ON OXYGEN-RELEASE BEHAVIOR\n');
fprintf('====================================================\n');

release_results = [];

for i = 1:length(x_values)
    for j = 1:length(y_values)
        x_current = x_values(i);
        y_current = y_values(j);
        M_formula_current = (1-x_current)*M_La + x_current*M_Ca + (1-y_current)*M_Co + y_current*M_Fe + 3*M_O;
        Delta_delta_current = 0.10 + 0.10*x_current + 0.08*y_current;
        Delta_delta_current = min(Delta_delta_current,0.20);
        eta_current = 0.90 - 0.05*x_current - 0.03*y_current;
        eta_current = max(0,min(eta_current,1));
        O2_release_wtpercent_current = (Delta_delta_current*M_O/M_formula_current)*100;
        OEC_cm3_g_current = ((eta_current*Delta_delta_current)/2)/M_formula_current*22414;
        release_results = [release_results; x_current, y_current, M_formula_current, Delta_delta_current, eta_current, O2_release_wtpercent_current, OEC_cm3_g_current];
    end
end

[max_release_wtpercent,idx_max_release] = max(release_results(:,6));
release_best_x = release_results(idx_max_release,1);
release_best_y = release_results(idx_max_release,2);

fprintf('Highest oxygen-release composition:\n');
fprintf('x = %.2f, y = %.2f\n',release_best_x,release_best_y);
fprintf('O2 release = %.4f wt%%\n',max_release_wtpercent);
fprintf('OEC        = %.3f cm^3(STP)/g\n',release_results(idx_max_release,7));

x_tick_step_rel = max(x_values)/7;                 % ~8 tick marks regardless of grid resolution
x_ticks_sparse_rel = 0:x_tick_step_rel:max(x_values);

fig = figure('Name','Doping Influence on Oxygen Release');
release_colors = lines(length(y_values));
hold on;
for j = 1:length(y_values)
    row_idx = release_results(:,2) == y_values(j);
    plot(release_results(row_idx,1),release_results(row_idx,6),'-o','LineWidth',1.8,'MarkerSize',4,'Color',release_colors(j,:),'MarkerFaceColor',release_colors(j,:),'DisplayName',sprintf('y = %.2f',y_values(j)));
end
hold off;
xlabel('Ca A-site doping, x','FontName','Times New Roman','FontSize',20,'FontWeight','bold');
ylabel('Oxygen release (wt%)','FontName','Times New Roman','FontSize',20,'FontWeight','bold');
title('Effect of A/B-Site Doping on Oxygen-Release Behavior','FontName','Times New Roman','FontSize',18,'FontWeight','bold');
xticks(x_ticks_sparse_rel);
legend('show','Location','eastoutside','NumColumns',2,'FontName','Times New Roman','FontSize',9,'FontWeight','bold');
ax = gca; ax.FontName = 'Times New Roman'; ax.FontSize = 16; ax.FontWeight = 'bold';
drawnow;
exportgraphics(fig,'Doping_Influence_Oxygen_Release.png','Resolution',1000);

Oxygen_Release_TPD = table(T_release_onset,T_release_peak,T_release_end,total_O2_released_mol_per_g,total_O2_released_cm3_g,total_O2_released_wtpercent,'VariableNames',{'Onset_Temp_C','Peak_Temp_C','End_Temp_C','O2_Released_mol_g','O2_Released_cm3STP_g','O2_Released_wtpercent'});

Oxygen_Reoxidation = table(delta_0,delta_max,eta_reversibility,reversible_delta,irreversible_delta,delta_reoxidized_final,irreversible_delta_loss,'VariableNames',{'Delta0','Delta_Max','Reversibility_Fraction','Reversible_Delta_Delta','Irreversible_Delta_Delta','Delta_After_Reox','Irreversible_Delta_Loss'});

Oxygen_Exchange_Capacity = table(reversible_delta,irreversible_delta,n_O2_exchange_mol_g,OEC_mmol_g,OEC_cm3_g,OEC_wtpercent,irreversible_wtpercent,'VariableNames',{'Reversible_Delta_Delta','Irreversible_Delta_Delta','OEC_mol_O2_g','OEC_mmol_O2_g','OEC_cm3STP_g','OEC_wtpercent','Irreversible_wtpercent'});

Doping_Oxygen_Release = array2table(release_results,'VariableNames',{'Ca_x','Fe_y','Formula_Molar_Mass_g_mol','Delta_Delta','Reversibility_Fraction','O2_Release_wtpercent','OEC_cm3STP_g'});


%% ============================================================
% STEP 19 — REDOX REACTION ENTHALPY (BEST FORMULATION)
% Redox reaction -> Reaction enthalpy
%
% NOTE: Delta_H_red is an ASSUMED linear doping model, not a fitted
% or measured value. Real values require isothermal TGA/DSC redox
% cycling or DFT formation-energy calculations. Typical literature
% ranges for Co/Fe-based perovskite oxides span roughly
% 120-400 kJ/mol O2, which is why the model is clamped to that range.
%% ============================================================

fprintf('\n====================================================\n');
fprintf('       REDOX REACTION ENTHALPY (SIMULATED)\n');
fprintf('====================================================\n');

Hred_base = 320.0;   % kJ/mol O2, assumed parent (undoped) LaCoO3-type enthalpy
Hred_best = Hred_base - 180*best_x - 120*best_y;
Hred_best = max(120.0,min(Hred_best,400.0));

fprintf('Assumed parent enthalpy (x=0,y=0) = %.1f kJ/mol O2\n',Hred_base);
fprintf('Assumed Delta_H_red (best formulation) = %.1f kJ/mol O2\n',Hred_best);
fprintf('NOTE: linear doping model, clamped to [120, 400] kJ/mol O2; not a measured value.\n');
fprintf('Replace with experimental DSC/TGA redox enthalpy or DFT formation energy before publication.\n');


%% ============================================================
% STEP 20 — REVERSIBLE OXYGEN EXCHANGE CAPACITY (STORAGE BASIS)
% Reaction enthalpy -> Oxygen exchange
%
% Reuses reversible_delta / n_O2_exchange_mol_g already computed in
% Step 17 (from eta_reversibility * Delta_delta_best) as the input
% to the storage-density calculation below.
%% ============================================================

fprintf('\n====================================================\n');
fprintf('   REVERSIBLE OXYGEN EXCHANGE CAPACITY (STORAGE BASIS)\n');
fprintf('====================================================\n');
fprintf('Reversible Delta-delta   = %.5f\n',reversible_delta);
fprintf('n(O2 exchanged)          = %.5e mol O2/g\n',n_O2_exchange_mol_g);
fprintf('Predicted OEC            = %.4f mmol O2/g = %.3f cm^3(STP)/g = %.4f wt%%\n',OEC_mmol_g,OEC_cm3_g,OEC_wtpercent);


%% ============================================================
% STEP 21 — THERMOCHEMICAL ENERGY-STORAGE DENSITY
% Oxygen exchange -> Storage capacity
%
% Energy_density = n(O2 exchanged, mol/kg) x Delta_H_red (kJ/mol O2)
%
% NOTE: This is a SIMULATED/PREDICTED value, not an experimentally
% demonstrated storage capacity — it is arithmetically consistent
% given the simulated OEC and the assumed (not measured) Delta_H_red
% from Step 19. Report as "predicted" until Delta_H_red and
% reversible_delta are replaced with measured DSC/TGA or DFT values.
%% ============================================================

fprintf('\n====================================================\n');
fprintf('       THERMOCHEMICAL ENERGY-STORAGE DENSITY\n');
fprintf('====================================================\n');

EnergyDensity_kJ_kg_best = n_O2_exchange_mol_g*1000*Hred_best;   % mol/kg x kJ/mol O2
EnergyDensity_Wh_kg_best = EnergyDensity_kJ_kg_best/3.6;
EnergyDensity_irrev_kJ_kg_best = ((irreversible_delta/2)/M_formula)*1000*Hred_best;

fprintf('Predicted reversible energy-storage density = %.3f kJ/kg  (= %.3f Wh/kg)\n',EnergyDensity_kJ_kg_best,EnergyDensity_Wh_kg_best);
fprintf('Predicted irreversible (lost) energy        = %.3f kJ/kg\n',EnergyDensity_irrev_kJ_kg_best);
fprintf('NOTE: based on the SIMULATED Delta_H_red and reversible_delta above; not an experimentally\n');
fprintf('demonstrated storage capacity. Report as predicted/simulated until validated against\n');
fprintf('measured DSC/TGA redox enthalpy or DFT formation energies.\n');


%% ============================================================
% STEP 22 — OPERATING TEMPERATURE RANGE
% Estimate the practical charging (reduction) and discharging
% (re-oxidation) temperature window.
%
% NOTE: The model has no kinetic hysteresis term, so the
% re-oxidation window is reported as the SAME temperature span as
% reduction (T_release_onset..T_release_end from Step 15), traversed
% in reverse on cooling. Real materials typically show a lower
% re-oxidation onset (kinetic hysteresis) that this simulation does
% not capture.
%% ============================================================

fprintf('\n====================================================\n');
fprintf('              OPERATING TEMPERATURE RANGE\n');
fprintf('====================================================\n');

fprintf('Charging (reduction) window   : %.1f - %.1f degC (peak %.1f degC)\n',T_release_onset,T_release_end,T_release_peak);
fprintf('Discharging (re-oxidation) window (idealized, no kinetic hysteresis): %.1f - %.1f degC\n',T_release_onset,T_release_end);
fprintf('Total operating span          = %.1f degC\n',T_release_end-T_release_onset);

fig = figure('Name','Operating Temperature Window');
barh(1,[T_release_onset, T_release_peak-T_release_onset, T_release_end-T_release_peak],'stacked');
hold on;
barh(2,[T_release_onset, T_release_peak-T_release_onset, T_release_end-T_release_peak],'stacked');
hold off;
colormap(gca,[1 1 1; 0.85 0.10 0.10; 0.93 0.49 0.19]);
set(gca,'YTick',[1 2],'YTickLabel',{'Charging (reduction)','Discharging (re-oxidation)'});
xlabel('Temperature (°C)','FontName','Times New Roman','FontSize',20,'FontWeight','bold');
title('Operating Temperature Window','FontName','Times New Roman','FontSize',20,'FontWeight','bold');
ax = gca; ax.FontName = 'Times New Roman'; ax.FontSize = 20; ax.FontWeight = 'bold';
drawnow;
exportgraphics(fig,'Operating_Temperature_Window.png','Resolution',1000);


%% ============================================================
% STEP 23 — CYCLING STABILITY MODEL (SIMULATED)
%
% NOTE: capacity_retention(N) = 100*exp(-k*N) with an ASSUMED,
% doping-dependent decay rate k. This is a placeholder for real
% redox-cycling TGA data, not a fitted or measured degradation
% curve. All retention numbers below are PREDICTED, not experimental.
%% ============================================================

fprintf('\n====================================================\n');
fprintf('           CYCLING STABILITY (SIMULATED)\n');
fprintf('====================================================\n');

decayrate_best = 0.05 + 0.02*best_x + 0.01*best_y;   % percent capacity loss per cycle (assumed)
cycle_numbers = (0:5:200)';
capacity_retention_best = 100*exp(-(decayrate_best/100)*cycle_numbers);

fprintf('Assumed decay rate (best formulation) = %.4f %%/cycle\n',decayrate_best);
fprintf('Predicted capacity retention after  50 cycles   = %.3f %%\n',100*exp(-(decayrate_best/100)*50));
fprintf('Predicted capacity retention after 100 cycles   = %.3f %%\n',100*exp(-(decayrate_best/100)*100));
fprintf('Predicted capacity retention after 200 cycles   = %.3f %%\n',100*exp(-(decayrate_best/100)*200));
fprintf('NOTE: replace this exponential decay model with real multi-cycle TGA/redox retention data.\n');

% Reference compositions for comparison: undoped parent vs best formulation
decayrate_parent = 0.05 + 0.02*0 + 0.01*0;
capacity_retention_parent = 100*exp(-(decayrate_parent/100)*cycle_numbers);

fig = figure('Name','Cycling Stability');
plot(cycle_numbers,capacity_retention_parent,'--','LineWidth',2,'Color',[0.50 0.50 0.50]);
hold on;
plot(cycle_numbers,capacity_retention_best,'LineWidth',2,'Color',[0.85 0.10 0.10]);
hold off;
xlabel('Cycle number','FontName','Times New Roman','FontSize',20,'FontWeight','bold');
ylabel('Predicted capacity retention (%)','FontName','Times New Roman','FontSize',20,'FontWeight','bold');
title('Cycling Stability','FontName','Times New Roman','FontSize',20,'FontWeight','bold');
legend({'Parent (x=0,y=0)','Best formulation'},'FontName','Times New Roman','FontSize',14,'FontWeight','bold');
ax = gca; ax.FontName = 'Times New Roman'; ax.FontSize = 20; ax.FontWeight = 'bold';
drawnow;
exportgraphics(fig,'Cycling_Stability.png','Resolution',1000);


%% ============================================================
% STEP 24 — STRUCTURE-VACANCY-STORAGE CORRELATION GRID
% Doping -> Structure -> Oxygen vacancies -> O2 exchange -> Storage
%
% Builds one row per (x,y) doping composition with structural
% (crystallite size), vacancy, exchange, storage-density and
% cycling-stability metrics, all under the SAME simulated models
% used above, so they can be correlated against each other.
%% ============================================================

K_scherrer = 0.9;
lambda_Cu = 0.15406;                 % nm
theta_ref = deg2rad(peak_positions(1)/2);

SVS_results = [];

for i = 1:length(x_values)
    for j = 1:length(y_values)
        x_c = x_values(i);
        y_c = y_values(j);

        M_formula_c = (1-x_c)*M_La + x_c*M_Ca + (1-y_c)*M_Co + y_c*M_Fe + 3*M_O;

        % Structural characteristic: XRD peak broadens with doping
        % (assumed lattice-strain proxy) -> smaller crystallite size
        peak_width_c = peak_width*(1 + 0.30*x_c + 0.20*y_c);
        crystallite_c = K_scherrer*lambda_Cu/(deg2rad(peak_width_c)*cos(theta_ref));

        Delta_delta_c = 0.10 + 0.10*x_c + 0.08*y_c;
        Delta_delta_c = min(Delta_delta_c,0.20);
        vacancy_percent_c = (delta_0+Delta_delta_c)/3*100;

        eta_c = 0.90 - 0.05*x_c - 0.03*y_c;
        eta_c = max(0,min(eta_c,1));
        reversible_delta_c = eta_c*Delta_delta_c;

        O2_release_wtpercent_c = (Delta_delta_c*M_O/M_formula_c)*100;
        OEC_wtpercent_c = (reversible_delta_c*M_O/M_formula_c)*100;

        Hred_c = Hred_base - 180*x_c - 120*y_c;
        Hred_c = max(120.0,min(Hred_c,400.0));

        n_O2_c = (reversible_delta_c/2)/M_formula_c;
        EnergyDensity_kJkg_c = n_O2_c*1000*Hred_c;
        EnergyDensity_Whkg_c = EnergyDensity_kJkg_c/3.6;

        decayrate_c = 0.05 + 0.02*x_c + 0.01*y_c;
        retention100_c = 100*exp(-(decayrate_c/100)*100);

        SVS_results = [SVS_results; x_c, y_c, M_formula_c, crystallite_c, vacancy_percent_c, O2_release_wtpercent_c, OEC_wtpercent_c, Hred_c, EnergyDensity_kJkg_c, EnergyDensity_Whkg_c, decayrate_c, retention100_c];
    end
end

fprintf('\n====================================================\n');
fprintf('   STRUCTURE-VACANCY-STORAGE CORRELATION (GRID)\n');
fprintf('====================================================\n');
fprintf('%-6s %-6s %-12s %-10s %-10s %-9s %-10s %-9s\n','x','y','Cryst(nm)','Vo(%)','O2rel(%)','OEC(%)','Ed(Wh/kg)','Ret100(%)');
for r = 1:size(SVS_results,1)
    fprintf('%-6.2f %-6.2f %-12.3f %-10.3f %-10.4f %-9.4f %-10.3f %-9.3f\n', ...
        SVS_results(r,1),SVS_results(r,2),SVS_results(r,4),SVS_results(r,5),SVS_results(r,6),SVS_results(r,7),SVS_results(r,10),SVS_results(r,12));
end


%% ============================================================
% STEP 25 — CORRELATION COEFFICIENTS
% Correlate structural characteristics with oxygen-vacancy
% formation, and relate vacancies to O2 exchange and storage.
%% ============================================================

r_cryst_vacancy = corrcoef(SVS_results(:,4),SVS_results(:,5));
r_vacancy_release = corrcoef(SVS_results(:,5),SVS_results(:,6));
r_vacancy_OEC = corrcoef(SVS_results(:,5),SVS_results(:,7));
r_OEC_energy = corrcoef(SVS_results(:,7),SVS_results(:,10));
r_energy_retention = corrcoef(SVS_results(:,10),SVS_results(:,12));

fprintf('\n====================================================\n');
fprintf('              CORRELATION COEFFICIENTS\n');
fprintf('====================================================\n');
fprintf('Crystallite size  vs  Oxygen vacancy (%%)   : r = %.4f\n',r_cryst_vacancy(1,2));
fprintf('Oxygen vacancy (%%) vs  O2 release (wt%%)     : r = %.4f\n',r_vacancy_release(1,2));
fprintf('Oxygen vacancy (%%) vs  OEC (wt%%)            : r = %.4f\n',r_vacancy_OEC(1,2));
fprintf('OEC (wt%%)          vs  Energy density (Wh/kg): r = %.4f\n',r_OEC_energy(1,2));
fprintf('Energy density      vs  Cycling retention    : r = %.4f\n',r_energy_retention(1,2));
fprintf('NOTE: correlations are computed across the SIMULATED doping grid, not experimental data.\n');

fig = figure('Name','Structure-Vacancy Correlation');
scatter(SVS_results(:,4),SVS_results(:,5),80,'filled','MarkerFaceColor',[0.00 0.45 0.74]);
xlabel('Crystallite size (nm)','FontName','Times New Roman','FontSize',20,'FontWeight','bold');
ylabel('Oxygen vacancy (%)','FontName','Times New Roman','FontSize',20,'FontWeight','bold');
title(sprintf('Crystallite Size vs Oxygen Vacancy (r = %.3f)',r_cryst_vacancy(1,2)),'FontName','Times New Roman','FontSize',20,'FontWeight','bold');
ax = gca; ax.FontName = 'Times New Roman'; ax.FontSize = 20; ax.FontWeight = 'bold';
drawnow;
exportgraphics(fig,'Structure_Vacancy_Correlation.png','Resolution',1000);

fig = figure('Name','Vacancy-OEC-EnergyDensity Correlation');
subplot(1,2,1);
scatter(SVS_results(:,5),SVS_results(:,7),70,'filled','MarkerFaceColor',[0.85 0.33 0.10]);
xlabel('Oxygen vacancy (%)','FontName','Times New Roman','FontSize',14,'FontWeight','bold');
ylabel('OEC (wt%)','FontName','Times New Roman','FontSize',14,'FontWeight','bold');
title(sprintf('r = %.3f',r_vacancy_OEC(1,2)),'FontName','Times New Roman','FontSize',14,'FontWeight','bold');
ax = gca; ax.FontName = 'Times New Roman'; ax.FontSize = 14; ax.FontWeight = 'bold';
subplot(1,2,2);
scatter(SVS_results(:,7),SVS_results(:,10),70,'filled','MarkerFaceColor',[0.47 0.67 0.19]);
xlabel('OEC (wt%)','FontName','Times New Roman','FontSize',14,'FontWeight','bold');
ylabel('Energy density (Wh/kg)','FontName','Times New Roman','FontSize',14,'FontWeight','bold');
title(sprintf('r = %.3f',r_OEC_energy(1,2)),'FontName','Times New Roman','FontSize',14,'FontWeight','bold');
ax = gca; ax.FontName = 'Times New Roman'; ax.FontSize = 14; ax.FontWeight = 'bold';
sgtitle('Oxygen Vacancy -> OEC -> Energy Density','FontName','Times New Roman','FontSize',18,'FontWeight','bold');
drawnow;
exportgraphics(fig,'Vacancy_OEC_EnergyDensity_Correlation.png','Resolution',1000);


%% ============================================================
% STEP 26 — THERMOCHEMICAL STORAGE PERFORMANCE COMPARISON

energy_matrix_grid = reshape(SVS_results(:,10),length(y_values),length(x_values));
retention_matrix_grid = reshape(SVS_results(:,12),length(y_values),length(x_values));
line_colors_storage = lines(length(y_values));
x_tick_step = max(x_values)/7;                    % ~8 tick marks regardless of grid resolution
x_ticks_sparse = 0:x_tick_step:max(x_values);

fig = figure('Name','Storage Performance - Energy Density');
hold on;
for j = 1:length(y_values)
    plot(x_values,energy_matrix_grid(j,:),'-o','LineWidth',1.8,'MarkerSize',4, ...
        'Color',line_colors_storage(j,:),'MarkerFaceColor',line_colors_storage(j,:), ...
        'DisplayName',sprintf('y = %.2f',y_values(j)));
end
hold off;
xlabel('Ca A-site doping, x','FontName','Times New Roman','FontSize',20,'FontWeight','bold');
ylabel('Predicted energy density (Wh/kg)','FontName','Times New Roman','FontSize',20,'FontWeight','bold');
title('Thermochemical Storage: Predicted Energy Density vs Doping','FontName','Times New Roman','FontSize',18,'FontWeight','bold');
xticks(x_ticks_sparse);
legend('show','Location','eastoutside','NumColumns',2,'FontName','Times New Roman','FontSize',9,'FontWeight','bold');
ax = gca; ax.FontName = 'Times New Roman'; ax.FontSize = 16; ax.FontWeight = 'bold';
drawnow;
exportgraphics(fig,'Storage_Performance_EnergyDensity.png','Resolution',1000);

fig = figure('Name','Storage Performance - Cycling Retention');
hold on;
for j = 1:length(y_values)
    plot(x_values,retention_matrix_grid(j,:),'-o','LineWidth',1.8,'MarkerSize',4, ...
        'Color',line_colors_storage(j,:),'MarkerFaceColor',line_colors_storage(j,:), ...
        'DisplayName',sprintf('y = %.2f',y_values(j)));
end
hold off;
xlabel('Ca A-site doping, x','FontName','Times New Roman','FontSize',20,'FontWeight','bold');
ylabel('Predicted retention after 100 cycles (%)','FontName','Times New Roman','FontSize',20,'FontWeight','bold');
title('Thermochemical Storage: Predicted Cycling Retention vs Doping','FontName','Times New Roman','FontSize',18,'FontWeight','bold');
xticks(x_ticks_sparse);
legend('show','Location','eastoutside','NumColumns',2,'FontName','Times New Roman','FontSize',9,'FontWeight','bold');
ax = gca; ax.FontName = 'Times New Roman'; ax.FontSize = 16; ax.FontWeight = 'bold';
drawnow;
exportgraphics(fig,'Storage_Performance_Retention.png','Resolution',1000);


%% ============================================================
% STEP 27 — OPTIMUM A/B-SITE DOPING COMPOSITION
% (STRUCTURE-VACANCY-STORAGE COMBINED SCORE)
%
% Combines vacancy formation, O2 exchange capacity, energy-storage
% density and cycling-stability retention (equal weights) into one
% normalized score across the SVS_results grid, independent of the
% simpler four-objective score computed in Step 1.
%% ============================================================

vacancy_col = SVS_results(:,5);
OEC_col = SVS_results(:,7);
energy_col = SVS_results(:,10);
retention_col = SVS_results(:,12);

vacancy_n = (vacancy_col-min(vacancy_col))/(max(vacancy_col)-min(vacancy_col));
OEC_n = (OEC_col-min(OEC_col))/(max(OEC_col)-min(OEC_col));
energy_n = (energy_col-min(energy_col))/(max(energy_col)-min(energy_col));
retention_n = (retention_col-min(retention_col))/(max(retention_col)-min(retention_col));

w_vac2 = 0.25; w_oec2 = 0.25; w_energy2 = 0.25; w_retention2 = 0.25;
combined_score = w_vac2*vacancy_n + w_oec2*OEC_n + w_energy2*energy_n + w_retention2*retention_n;

[best_combined_score,idx_combined] = max(combined_score);
combined_best_x = SVS_results(idx_combined,1);
combined_best_y = SVS_results(idx_combined,2);

fprintf('\n====================================================\n');
fprintf('   OPTIMUM DOPING - STRUCTURE-VACANCY-STORAGE SCORE\n');
fprintf('====================================================\n');
fprintf('Weights: vacancy=%.2f, OEC=%.2f, energy=%.2f, retention=%.2f\n',w_vac2,w_oec2,w_energy2,w_retention2);
fprintf('Grid tested: x = %.2f:%.2f:%.2f, y = %.2f:%.2f:%.2f\n',min(x_values),x_values(2)-x_values(1),max(x_values),min(y_values),y_values(2)-y_values(1),max(y_values));
fprintf('Optimum x = %.2f\n',combined_best_x);
fprintf('Optimum y = %.2f\n',combined_best_y);
fprintf('Combined score (0-1 normalized) = %.4f\n',best_combined_score);
fprintf('Crystallite size     = %.3f nm\n',SVS_results(idx_combined,4));
fprintf('Oxygen vacancy       = %.4f %%\n',SVS_results(idx_combined,5));
fprintf('O2 release           = %.4f wt%%\n',SVS_results(idx_combined,6));
fprintf('OEC                  = %.4f wt%%\n',SVS_results(idx_combined,7));
fprintf('Predicted energy density       = %.3f Wh/kg\n',SVS_results(idx_combined,10));
fprintf('Predicted retention (100 cyc)  = %.3f %%\n',SVS_results(idx_combined,12));
if combined_best_x == best_x && combined_best_y == best_y
    fprintf('\nConsistent with the Step 1 optimization result (x=%.2f, y=%.2f).\n',best_x,best_y);
else
    fprintf('\nNOTE: differs from the Step 1 optimization result (x=%.2f, y=%.2f) because\n',best_x,best_y);
    fprintf('this score additionally weights energy density and cycling stability.\n');
end
if combined_best_x == max(x_values) || combined_best_y == max(y_values)
    fprintf('\nNOTE: optimum touches the upper boundary of the tested x/y range (%.2f) on at\n',max(x_values));
    fprintf('least one axis; extend the grid further to confirm this is a true optimum.\n');
end

fprintf('\n====================================================\n');
fprintf('        OPTIMIZATION CRITERIA COMPARISON (A vs B)\n');
fprintf('====================================================\n');
fprintf('Score A (Step 1, "optimization_score"):\n');
fprintf('  Inputs : vacancy proxy, OEC proxy, ENERGY-DENSITY PROXY (a simple linear\n');
fprintf('           function of x,y — NOT the physically derived storage density),\n');
fprintf('           and a stability proxy.\n');
fprintf('  Weights: vacancy=%.2f, OEC=%.2f, energy=%.2f, stability=%.2f\n',w_vacancy,w_oec,w_energy,w_stability);
fprintf('  Result : x = %.2f, y = %.2f (score = %.4f)\n',best_x,best_y,best_score);
fprintf('\nScore B (Step 27, "combined_score", Structure-Vacancy-Storage):\n');
fprintf('  Inputs : vacancy percent, OEC (wt%%), the PHYSICALLY DERIVED predicted\n');
fprintf('           energy density (mol O2 exchanged x assumed Delta_H_red), and\n');
fprintf('           predicted cycling retention after 100 cycles. No stability term.\n');
fprintf('  Weights: vacancy=%.2f, OEC=%.2f, energy=%.2f, retention=%.2f\n',w_vac2,w_oec2,w_energy2,w_retention2);
fprintf('  Result : x = %.2f, y = %.2f (score = %.4f)\n',combined_best_x,combined_best_y,best_combined_score);
fprintf('\nThese are DIFFERENT optimization criteria (different objective sets and\n');
fprintf('weights), so they are not required to agree, and here they do not: Score A\n');
fprintf('optimizes for x=%.2f,y=%.2f while Score B optimizes for x=%.2f,y=%.2f.\n',best_x,best_y,combined_best_x,combined_best_y);
fprintf('Score A''s "energy" term is a simple linear proxy used only for the initial\n');
fprintf('doping-matrix screen; Score B''s energy term is the physically derived\n');
fprintf('predicted storage density and additionally accounts for cycling stability,\n');
fprintf('which Score A does not consider at all. For reporting, present both optima\n');
fprintf('explicitly labeled by criterion (Score A vs Score B) rather than a single\n');
fprintf('"the" optimum — e.g.: "Under the screening criterion (Score A), the optimum\n');
fprintf('composition is x=%.2f, y=%.2f; under the structure-vacancy-storage criterion\n',best_x,best_y);
fprintf('(Score B), which additionally weights predicted energy density and cycling\n');
fprintf('retention, the optimum is x=%.2f, y=%.2f."\n',combined_best_x,combined_best_y);

Optimization_Criteria_Comparison = table(best_x,best_y,best_score,combined_best_x,combined_best_y,best_combined_score, ...
    'VariableNames',{'ScoreA_Best_x','ScoreA_Best_y','ScoreA_Value','ScoreB_Best_x','ScoreB_Best_y','ScoreB_Value'});

Reaction_Enthalpy_Storage = table(Hred_base,Hred_best,EnergyDensity_kJ_kg_best,EnergyDensity_Wh_kg_best,EnergyDensity_irrev_kJ_kg_best,T_release_onset,T_release_peak,T_release_end,'VariableNames',{'Hred_Parent_kJ_mol_O2','Hred_Best_kJ_mol_O2','EnergyDensity_kJ_kg_Predicted','EnergyDensity_Wh_kg_Predicted','Irreversible_Energy_kJ_kg_Predicted','Charging_Onset_C','Charging_Peak_C','Charging_End_C'});

Cycling_Stability = table(cycle_numbers,capacity_retention_parent,capacity_retention_best,'VariableNames',{'Cycle_Number','Retention_Parent_Percent_Predicted','Retention_Best_Percent_Predicted'});

Structure_Vacancy_Storage = array2table(SVS_results,'VariableNames',{'Ca_x','Fe_y','Formula_Molar_Mass_g_mol','Crystallite_Size_nm','Vacancy_Percent','O2_Release_wtpercent','OEC_wtpercent','Hred_kJ_mol_O2','EnergyDensity_kJ_kg_Predicted','EnergyDensity_Wh_kg_Predicted','Decay_Rate_Percent_per_Cycle','Retention_100cyc_Percent_Predicted'});
Structure_Vacancy_Storage.Combined_Score = combined_score;

Correlation_Coefficients = table(r_cryst_vacancy(1,2),r_vacancy_release(1,2),r_vacancy_OEC(1,2),r_OEC_energy(1,2),r_energy_retention(1,2),'VariableNames',{'r_Crystallite_Vacancy','r_Vacancy_O2Release','r_Vacancy_OEC','r_OEC_EnergyDensity','r_EnergyDensity_Retention'});


%% ============================================================
% STEP 28 — SAVE THE IMPORTANT MATRICES TO ONE EXCEL FILE
%
% Only the key doping-grid MATRICES plus a one-row best-formulation
% summary are written here — the old version dumped 30 sheets
% including raw spectral curves (XRD/BET traces, SEM particle list,
% TGA curve, etc.). Those raw-curve tables are no longer written.
% Everything below is still computed earlier in the script; this
% section only controls what actually gets saved to disk.
%% ============================================================

excel_file = 'Best_Formulation_Doping_Matrices.xlsx';

% --- (1) One-row summary of the best formulation (both scoring criteria) ---
Best_Formulation = table(best_x,best_y,best_score,string(best_formula),crystallite_size, ...
    BET_surface_area,total_pore_volume,average_pore_diameter,delta_0,Delta_delta_best, ...
    delta_absolute(end),vacancy_percentage(end),Hred_best,EnergyDensity_kJ_kg_best, ...
    EnergyDensity_Wh_kg_best,decayrate_best,combined_best_x,combined_best_y,best_combined_score, ...
    'VariableNames',{'ScoreA_Best_x','ScoreA_Best_y','ScoreA_Value','Formulation','Crystallite_Size_nm', ...
    'BET_Surface_Area_m2_g','Total_Pore_Volume_cm3_g','Average_Pore_Diameter_nm','Delta0_Simulated_Reference', ...
    'Maximum_Delta_Delta','Final_Absolute_Delta','Final_Oxygen_Site_Deficiency_Percent','Hred_kJ_mol_O2_Assumed', ...
    'EnergyDensity_kJ_kg_Predicted','EnergyDensity_Wh_kg_Predicted','Decay_Rate_Percent_per_Cycle_Assumed', ...
    'ScoreB_Best_x','ScoreB_Best_y','ScoreB_Value'});

% --- (2) Step 1 doping-matrix screening grid (Score A, normalized components) ---
Optimization_Matrix = table(X(:),Y(:),Vo_norm(:),OEC_norm(:),Energy_norm(:),Stability_norm(:),optimization_score(:), ...
    'VariableNames',{'Ca_x','Fe_y','Vo_Normalized','OEC_Normalized','EnergyDensity_Normalized','Stability_Normalized','ScoreA_Optimization_Score'});

% --- (3) A/B-site doping vs oxygen-vacancy grid ---
AB_Doping_Vacancy_Matrix = array2table(results_AB,'VariableNames',{'Ca_x','Fe_y','Delta_Delta','Absolute_Delta','Vacancy_Fraction','Vacancy_Percent'});

% --- (4) A/B-site doping vs oxygen-release grid ---
Doping_Oxygen_Release_Matrix = array2table(release_results,'VariableNames',{'Ca_x','Fe_y','Formula_Molar_Mass_g_mol','Delta_Delta','Reversibility_Fraction','O2_Release_wtpercent','OEC_cm3STP_g'});

% --- (5) Combined Structure-Vacancy-Storage grid (Score B) — the most complete matrix ---
Structure_Vacancy_Storage = array2table(SVS_results,'VariableNames',{'Ca_x','Fe_y','Formula_Molar_Mass_g_mol','Crystallite_Size_nm','Vacancy_Percent','O2_Release_wtpercent','OEC_wtpercent','Hred_kJ_mol_O2','EnergyDensity_kJ_kg_Predicted','EnergyDensity_Wh_kg_Predicted','Decay_Rate_Percent_per_Cycle','Retention_100cyc_Percent_Predicted'});
Structure_Vacancy_Storage.Combined_Score = combined_score;

% --- (6) Correlation coefficients between the matrices above ---
Correlation_Coefficients = table(r_cryst_vacancy(1,2),r_vacancy_release(1,2),r_vacancy_OEC(1,2),r_OEC_energy(1,2),r_energy_retention(1,2), ...
    'VariableNames',{'r_Crystallite_Vacancy','r_Vacancy_O2Release','r_Vacancy_OEC','r_OEC_EnergyDensity','r_EnergyDensity_Retention'});

% --- (7) Score A vs Score B optimum comparison ---
Optimization_Criteria_Comparison = table(best_x,best_y,best_score,combined_best_x,combined_best_y,best_combined_score, ...
    'VariableNames',{'ScoreA_Best_x','ScoreA_Best_y','ScoreA_Value','ScoreB_Best_x','ScoreB_Best_y','ScoreB_Value'});

% --- (8) Cycling-stability curve (parent vs best formulation) ---
Cycling_Stability_Matrix = table(cycle_numbers,capacity_retention_parent,capacity_retention_best, ...
    'VariableNames',{'Cycle_Number','Retention_Parent_Percent_Predicted','Retention_Best_Percent_Predicted'});

%% ============================================================
% WRITE EXCEL FILE (single workbook, important matrices only)
%% ============================================================

writetable(Best_Formulation,excel_file,'Sheet','Best_Formulation');
writetable(Optimization_Matrix,excel_file,'Sheet','Optimization_Matrix_ScoreA');
writetable(AB_Doping_Vacancy_Matrix,excel_file,'Sheet','AB_Doping_Vacancy_Matrix');
writetable(Doping_Oxygen_Release_Matrix,excel_file,'Sheet','Doping_O2Release_Matrix');
writetable(Structure_Vacancy_Storage,excel_file,'Sheet','Structure_Vacancy_Storage');
writetable(Correlation_Coefficients,excel_file,'Sheet','Correlation_Coefficients');
writetable(Optimization_Criteria_Comparison,excel_file,'Sheet','ScoreA_vs_ScoreB');
writetable(Cycling_Stability_Matrix,excel_file,'Sheet','Cycling_Stability');


%% ============================================================
% FINAL OUTPUT
%% ============================================================

fprintf('\n====================================================\n');
fprintf('             FINAL ANALYSIS SUMMARY\n');
fprintf('====================================================\n');
fprintf('\nBest formulation:\n');
fprintf('La_{%.2f}Ca_{%.2f}Co_{%.2f}Fe_{%.2f}O3\n',1-best_x,best_x,1-best_y,best_y);
fprintf('\nOxygen non-stoichiometry:\n');
fprintf('Simulated/Reference delta0 = %.5f\n',delta_0);
fprintf('Maximum Delta-delta     = %.5f\n',Delta_delta_best);
fprintf('Maximum absolute delta  = %.5f\n',delta_absolute(end));
fprintf('\nOxygen-site deficiency:\n');
fprintf('Vo/3 fraction             = %.5f\n',vacancy_fraction(end));
fprintf('Oxygen-site deficiency    = %.3f %% of the three nominal oxygen sites (delta = %.4f)\n',vacancy_percentage(end),delta_absolute(end));
fprintf('\nCo oxidation state:\n');
fprintf('Co2+ = %.2f %%\n',Co2_fraction*100);
fprintf('Co3+ = %.2f %%\n',Co3_fraction*100);
fprintf('Average Co valence = %.3f+\n',Co_average_valence);
fprintf('\nFe oxidation state:\n');
fprintf('Fe2+ = %.2f %%\n',Fe2_fraction*100);
fprintf('Fe3+ = %.2f %%\n',Fe3_fraction*100);
fprintf('Average Fe valence = %.3f+\n',Fe_average_valence);
fprintf('\nO 1s analysis:\n');
fprintf('Lattice O      = %.2f %%\n',lattice_oxygen_fraction*100);
fprintf('Vacancy-related O = %.2f %%\n',vacancy_related_fraction*100);
fprintf('Adsorbed O     = %.2f %%\n',adsorbed_oxygen_fraction*100);
fprintf('\nHighest vacancy formulation:\n');
fprintf('Ca x = %.2f\n',vacancy_x);
fprintf('Fe y = %.2f\n',vacancy_y);
fprintf('Absolute delta = %.5f\n',vacancy_delta);
fprintf('Vacancy = %.3f %%\n',maximum_vacancy);
fprintf('\nOxygen release (TPD):\n');
fprintf('Peak release temperature = %.1f degC\n',T_release_peak);
fprintf('Total O2 released        = %.4f wt%%\n',total_O2_released_wtpercent);
fprintf('\nReversible oxygen uptake:\n');
fprintf('Reversibility fraction   = %.3f\n',eta_reversibility);
fprintf('Delta after re-oxidation = %.5f\n',delta_reoxidized_final);
fprintf('\nOxygen exchange capacity (predicted):\n');
fprintf('OEC = %.4f wt%% = %.3f cm^3(STP)/g = %.4f mmol O2/g\n',OEC_wtpercent,OEC_cm3_g,OEC_mmol_g);
fprintf('\nHighest oxygen-release formulation:\n');
fprintf('Ca x = %.2f\n',release_best_x);
fprintf('Fe y = %.2f\n',release_best_y);
fprintf('O2 release = %.4f wt%%\n',max_release_wtpercent);
fprintf('\nReaction enthalpy / storage density (best formulation, ASSUMED/PREDICTED):\n');
fprintf('Assumed Delta_H_red = %.1f kJ/mol O2\n',Hred_best);
fprintf('Predicted energy-storage density = %.3f kJ/kg (%.3f Wh/kg)\n',EnergyDensity_kJ_kg_best,EnergyDensity_Wh_kg_best);
fprintf('Operating temperature range = %.1f - %.1f degC\n',T_release_onset,T_release_end);
fprintf('\nCycling stability (best formulation, PREDICTED):\n');
fprintf('Assumed decay rate = %.4f %%/cycle; predicted retention after 100 cycles = %.3f %%\n',decayrate_best,100*exp(-(decayrate_best/100)*100));
fprintf('\nTwo optimization criteria were evaluated (see OPTIMIZATION CRITERIA\n');
fprintf('COMPARISON section above for full explanation):\n');
fprintf('  Score A (Step 1 screening)              : x = %.2f, y = %.2f (score = %.4f)\n',best_x,best_y,best_score);
fprintf('  Score B (structure-vacancy-storage)      : x = %.2f, y = %.2f (score = %.4f)\n',combined_best_x,combined_best_y,best_combined_score);
if ~(combined_best_x==best_x && combined_best_y==best_y)
    fprintf('  These differ because Score B uses the physically derived predicted energy\n');
    fprintf('  density and predicted cycling retention instead of Score A''s linear proxies.\n');
end
if combined_best_x == max(x_values) || combined_best_y == max(y_values) || best_x == max(x_values) || best_y == max(y_values)
    fprintf('  NOTE: at least one optimum touches the grid boundary (%.2f) — extend the\n',max(x_values));
    fprintf('  grid further to confirm these are true optima.\n');
end
fprintf('\nExcel workbook saved as:\n');
fprintf('%s\n',excel_file);
fprintf('\nGenerated plots:\n');
fprintf('1. XRD_best_formulation.png\n');
fprintf('2. SEM_morphology.png\n');
fprintf('3. SEM_particle_size.png\n');
fprintf('4. BET_isotherm.png\n');
fprintf('5. Oxygen_NonStoichiometry.png\n');
fprintf('6. O1s_XPS_Oxygen_Species.png\n');
fprintf('7. Co_Fe_Oxidation_State.png\n');
fprintf('8. AB_Doping_Oxygen_Vacancy.png\n');
fprintf('9. Doping_vs_Oxygen_Vacancy.png\n');
fprintf('10. Oxygen_Release_TPD_Profile.png\n');
fprintf('11. Redox_Hysteresis_Release_Uptake.png\n');
fprintf('12. Oxygen_Exchange_Capacity.png\n');
fprintf('13. Doping_Influence_Oxygen_Release.png\n');
fprintf('14. Operating_Temperature_Window.png\n');
fprintf('15. Cycling_Stability.png\n');
fprintf('16. Structure_Vacancy_Correlation.png\n');
fprintf('17. Vacancy_OEC_EnergyDensity_Correlation.png\n');
fprintf('18. Storage_Performance_EnergyDensity.png\n');
fprintf('19. Storage_Performance_Retention.png\n');
fprintf('\n====================================================\n');
fprintf('             ANALYSIS COMPLETED\n');
fprintf('====================================================\n');
fprintf('\nIMPORTANT: This is a SIMULATION-STAGE workflow. XRD, SEM, BET,\n');
fprintf('TGA/TPD, XPS, delta0, and OEC values are all simulated/assumed\n');
fprintf('for workflow validation, not experimental results. Before\n');
fprintf('publication: (1) replace delta0 with a real iodometric protocol,\n');
fprintf('(2) fit real XPS spectra instead of assumed fractions, and\n');
fprintf('(3) validate the optimization against measured Vo, OEC, energy\n');
fprintf('density, and cycling stability rather than the assumed model.\n');
fprintf('(4) replace Delta_H_red with measured DSC/TGA redox enthalpy or DFT formation\n');
fprintf('energies, and (5) replace the assumed exponential cycling-decay model with\n');
fprintf('real multi-cycle TGA/redox retention data.\n');