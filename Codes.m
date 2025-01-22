%Basic lattice strain model- getting the prameters for the equation

clear all
clc


T=X+273; %Replace "X" with the temperature in celcius
r=[]; %radii of elements being modeled
known(2,2:length(r))=NaN;
known(1)=T;
known(2,:)=r;


PartCoef=[]; %put in known partition coefficients here
beta0 = [1 1 100];

[output,resid,J,sigma]=nlinfit(known,PartCoef,@lsmf,beta0); % @lsmf reads the function in the lsmf.m file. Ensure the lsmf.m and this codes.m file are saved in the same folder for this to work.
[Ypred,delta] = nlpredci(@lsmf,known,real(output),resid,'Covar',sigma,'alpha',0.33);
ci = nlparci(beta0,resid,'covar',sigma,'alpha',0.33);

%type "output" into the command window to get the modeled parameters for
%D0, r0, and E

%%
%Basic lattice strain model- getting predictions and plotting with modeled values

clear all
clc

D0=[];%modeled value from the previous analysis
r0=[];%modeled value from the previous analysis
E=[];%modeled value from the previous analysis
T=X+273;%Replace "X" with the temperature in celcius
r=[]; %radii of elements being modeled
unknown(1)=D0;
unknown(2)=r0;
unknown(3)=E;
known(2,2:length(r))=NaN;
known(1)=T;
known(2,:)=r;

PredictedD=lsmf(unknown,known);%this feeds in the modeled values into the model function and provides predicted partition coefficient values 

PartCoef=[]; %put in known partition coefficients here


%plot the predicted and actual data into a parabola
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
%'La' 'Ce' 'Nd' 'Sm' 'Eu' 'Gd' 'Tb' 'Dy' 'Er' 'Yb' 'Lu'
PXsource=[]'; %starting REE composition
DOl=[]'; %Olivine partition coefficients
DGt=[]';%Garnet partition coefficients
DOpx=[]'; %Orthopyroxene partition coefficients
DCpx=[]';%Clinopyroxene partition coefficients
%Chondrite values from Anders and Grevesse
AG=[0.2347 0.6032 0.4524 0.1471 0.1966 0.036 0.2427 0.1589 0.1625 0.0243]';
ElementOrder=(1:10);
%Melt fractions
F=[0.01 0.02 0.05 0.075 0.1 0.125 0.15 0.175 0.2 0.225 0.25 0.275 0.3];


%Instructions: for each mineral, include the % you want to assess (e.g.,
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

%Plot shows the change in elements over different melt fractions from the
%model
figure
semilogy(ElementOrder,LiquidBatch,'-','LineWidth',3);
hold on
axis square
title('Liquid Batch Melting of a Gt-bearing Pyroxenite');
ylabel('chondrite normalized');
set(gca,'TickLength',[0.03 0.03]);



%These models will give predicted concentrations for the trace elements
%analyzed. These predicted values can then be used for ratio comparisons (e.g., Ce/CM, Yb/Sm). 




