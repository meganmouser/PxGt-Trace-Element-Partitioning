%Lattice strain model fit: first define the input prameters for the equation

clear all
clc

%Example lattice strain model for REE3+ in garnet X-site

T=1175+273; %Replace "X" with the temperature in celcius
r=[1.16 1.143 1.126 1.109 1.079 1.053 1.04 1.027 1.015 1.004 0.994 0.985 0.977]; %radii of elements being modeled (in angstroms)
known(2,2:length(r))=NaN;
known(1)=T;
known(2,:)=r;


PartCoef=[0.0038 0.0060 0.0114 0.0238 0.0972 0.2853 0.5258 0.8671 1.3733 2.066 2.8755 3.7557 4.8705]; %put in known partition coefficients here
beta0 = [1 1 100];%initial parameter (D0, r0, E) predictions

[output,resid,J,sigma]=nlinfit(known,PartCoef,@lsmf,beta0); % @lsmf reads the function in the lsmf.m file. Ensure the lsmf.m and this codes.m file are saved in the same folder for this to work.
[Ypred,delta] = nlpredci(@lsmf,known,real(output),resid,'Covar',sigma,'alpha',0.33);
ci = nlparci(beta0,resid,'covar',sigma,'alpha',0.33);

%type "output" into the command window to get the modeled parameters for
%D0, r0, and E

%%
%Compare a lattice strain model fit to the experimental measurements

clear all
clc

D0=[54.8641];%modeled value from the previous analysis
r0=[0.8119];%modeled value from the previous analysis
E=[343.4178];%modeled value from the previous analysis
T=1175+273;%Replace "X" with the temperature in celcius
r=[1.16 1.143 1.126 1.109 1.079 1.053 1.04 1.027 1.015 1.004 0.994 0.985 0.977]; %radii of elements being modeled
unknown(1)=D0;
unknown(2)=r0;
unknown(3)=E;
known(2,2:length(r))=NaN;
known(1)=T;
known(2,:)=r;

PredictedD=lsmf(unknown,known);%this feeds in the modeled values into the model function and provides predicted partition coefficients 

PartCoef=[0.0038 0.0060 0.0114 0.0238 0.0972 0.2853 0.5258 0.8671 1.3733 2.066 2.8755 3.7557 4.8705]; %input known partition coefficients here


%plot the predicted and measured data 
figure
semilogy(r,PredictedD,'k-'); %line plot of predicted partition coefficients
hold on
plot(r,PartCoef,'ko');%known partition coefficients
axis square 
set(gca,'TickLength',[0.03;0.03]);
ylabel('Partition Coefficient');
xlabel('Radius (A)');
legend('Predicted','Measured');




%%
%Lunar Model

clear all
clc

%Example 

%'La' 'Ce' 'Nd' 'Sm' 'Gd' 'Tb' 'Dy' 'Er' 'Yb' 'Lu'
PXsource=[4.0815 11.8494 10.5714 3.9303 5.6605 1.0771 7.3218 4.8316 4.6878 0.6872]'; %starting REE composition
DOl=[0.0001 0.0001 0.0001 0.0006 0.001 0.002 0.003 0.008 0.019 0.03]'; %Olivine partition coefficients
DGt=[0.036 0.0091 0.0519 0.2166 0.6092 1.0083 1.5726 2.9704 4.4424 5.12]';%Garnet partition coefficients
DOpx=[0.007 0.009 0.014 0.022 0.037 0.048 0.06 0.1 0.17 0.22]'; %Orthopyroxene partition coefficients
DCpx=[0.02661 0.0437 0.08002 0.11696 0.14678 0.15146 0.15606 0.15815 0.16052 0.17007]';%Clinopyroxene partition coefficients 
%Chondrite values from Anders and Grevesse
AG=[0.2347 0.6032 0.4524 0.1471 0.1966 0.036 0.2427 0.1589 0.1625 0.0243]';
ElementOrder=(1:10);
%Melt fractions
F=[0.01 0.02 0.05 0.075 0.1 0.125 0.15 0.175 0.2 0.225 0.25 0.275 0.3];


%Instructions: for each mineral, include the degree of melting you want to assess (e.g.,
%10%=0.1), the sum of all % must be 1. If you wish to model without a given mineral, make sure to
%remove it from the "Bulk" equation or set it to zero. 
XDoliv=DOl.*0.2;
XDgt=DGt.*0.2;
XDcpx=DCpx.*0.2;
XDopx=DOpx.*0.4;
Bulk=XDoliv+XDcpx+XDgt+XDopx;

%modeling the inital composition, chondrite normalized
Cs0=PXsource./AG;
%modeling the liquid composition over a range of melt fractions (F)
LiquidBatch=Cs0./(Bulk+F.*(F-Bulk));

%Plot shows the change in elemental abundances over different melt fractions from the
%model
figure
semilogy(ElementOrder,LiquidBatch,'-','LineWidth',3);
hold on
axis square
title('Liquid Batch Melting of a Gt-bearing Pyroxenite');
ylabel('chondrite normalized');
%axis([1 10 10.^-1 10.^3]);
set(gca,'TickLength',[0.03 0.03]);
xticks([1 2 3 4 5 6 7 8 9 10]);
xticklabels({'La','Ce','Nd','Sm','Gd','Tb','Dy','Er','Yb','Lu'});



%This model provides predicted concentrations for the trace elements
%analyzed. These predicted values can then be used for ratio comparisons (e.g., Ce/Sm, Yb/Sm). 




