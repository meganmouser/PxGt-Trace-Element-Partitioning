function [D]=lsmf(unknown,known)
D0= unknown(1); 
r0= unknown(2).*10.^-10;
E= unknown(3).*10.^9;
T= known(1); %temperature, calls from the "code" file
r= known(2,:).*10.^-10; %known cation radii, calls from the "code" file
R=8.314; %gas constant
Na=6.022.*10.^23; %Avodagro's number
A=(-4.*pi.*E.*Na)./(R.*T);
B=(((r0./2).*(r0-r).^2)-(1/3.*(r0-r)).^3);
D= D0.*exp(A.*B);
end

