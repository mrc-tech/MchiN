close all
clear
clc

% addpath('..')

D  = 1.2; % diametro palo [m]
fc = 20e3; % resistenza calcestruzzo [kPa]
fy = 400e3; % resistenza acciaio [kPa]
db = 0.02; % diametro barre [m]
nb = 20; % numero barre
c  = 0.05; % copriferro [m]
rho = nb*(db/D)^2; % percentuale geometrica armatura
N = 10000; % sforzo normale [kN]

sez = MchiNstrisce3(D/1000); % definisce la sezione
sez.defineMaterials(fc, fy,  -0.002,-0.0035, 0.002,0.0675); % definisce i materiali
sez.initCirc(D, db, rho, c); % inizializza la sezione circolare

% [Mu, phiu, yn] = sez.ultimatePoint(-N); % compressione: negativa
[~, ~, Mu, phiu] = sez.findPoints2(-N, 'plotCurve');

figure
sez.plotSection
disp(['momento resistente : ' num2str(Mu)])
disp(['curvatura ultima   : ' num2str(phiu)])