% Starts a work session for the petide spectra alignment project

global AAletters AAmasses AAcodes AAnames AAsnps;
global bProfile yProfile bNames yNames;

load AAdata.mat;
load AAsnps;
load AAmasses_modC;   fprintf(1,'Warning: loaded AAmasses with modified Cysteine (+57 Da).\n');
load -mat msProfileOrbitrapETD.mat;
% load -mat msProfileOrbitrap.mat;
%load -mat msProfileReal.mat;
% load -mat msProfileISBinteger.mat;
%load -mat msProfileISB.mat;
%load -mat proteinData.dat;
%load -mat data_spectrin.mat;
