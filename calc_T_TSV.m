function T_TSV = calc_T_TSV(cond_si,esi,epoly,cond_cu,tsv_length,rvia,tox,pitch,freq)
%Copyright (c) 2018 Hakki Mert Torun
%Calculate S-Parameters of a coupled TSV array considering depletion
%capacitance. The S-Parameters are then converted to Y-Matrix, the matrix
%is reduced for assigning power/signal/ground, then 2 port T-Parameters are
%calculated.
%Impelementation based on the paper:
%Engin, A. Ege, and Srinidhi Raghavan Narasimhan.
%"Modeling of crosstalk in through silicon vias." 
%IEEE Trans. Electromagn. Compat. 55.1 (2013): 149-158.
%%
mu0 = 4*pi*1e-7;
e0 = 8.854*1e-12;
tdep = 0.8e-6;
%% Create a (rowNumber x columnNumber) grid for TSV array
columnNumber = 10;
rowNumber = 10;
%Index of Initial Ground TSV ([1,1] is top left corner)
indexOfRef = [1,1];
dn0 = zeros(rowNumber,columnNumber);

dn0(end:-1:1,:) = repmat(1j*(1:1:rowNumber)',1,columnNumber);
dn0(:,1:end) = dn0 + repmat(1:1:columnNumber,rowNumber,1);

dij = pdist2((1:1:rowNumber)',(1:1:columnNumber)');
dij = dij*(pitch);

dn0 = dn0-dn0(indexOfRef(1),indexOfRef(2));
dij_temp = reshape(dn0,1,rowNumber*columnNumber);

dij_temp(real(dij_temp) == 0 & imag(dij_temp) == 0) = [];

dij_r = pdist2(real(dij_temp)',real(dij_temp)').^2;
dij_im = pdist2(imag(dij_temp)',imag(dij_temp)').^2;

dij = sqrt(dij_r+dij_im);
dij = dij.*(pitch);
dij(dij == 0) = 1;
dn0 = abs(dn0)*(pitch);
dn0(dn0 == 0) = [];
%% Calculate Self Inductance of each element
Lii = mu0/2/pi*2*log(dn0/(rvia+tox+tdep));
L = dn0'*dn0;
L = mu0/2/pi*log(L/(rvia+tox+tdep));
L = L-mu0/2/pi*log(dij);
L(1:rowNumber*columnNumber:end) = Lii;

%%

%% Calculate Self & Mutual Capacitance and Inductance

C = mu0*e0*esi*inv(L);
% Cii = pi*e0/log(pitch/(rvia+tox+tdep));
% Cii = sum(C,2);
% C(logical(eye(size(C,1)))) = Cii;
G = mu0*cond_si*inv(L);
% Gii = sum(G,2);
% G(logical(eye(size(G,1)))) = Gii;

Lii = mu0/2/pi*2*log(dn0/(rvia));
L = dn0'*dn0;
L = mu0/2/pi*log(L/(rvia));
L = L-mu0/2/pi*log(dij);
L(1:rowNumber*columnNumber:end) = Lii;

%% Calculate Depletion Capacitance of each TSV
Cnn = (log((rvia+tox)/rvia)/(2*pi*epoly*e0)+log((rvia+tox+tdep)/(rvia+tox))/(2*pi*esi*e0))^-1;

% freq = linspace(0.1,20,200)*1e9;

w = 2*pi*freq;
Z = sqrt(1j*w*mu0/cond_cu)/(2*pi*(rvia)).*besseli(0,rvia*sqrt(1j*w*mu0*cond_cu))./besseli(1,rvia*sqrt(1j*w*mu0*cond_cu));
Z = 2*Z;
Lt = zeros(size(L,1),size(L,2),length(freq));
R = Lt;

Zt = Lt; 
Yt = Lt; 

%% Matrix Reduction
%Ground Elements: ground TSV assignment in (rowNumber x columnNumber) grid for TSV array
ground_elements = 2:2:(rowNumber*columnNumber-1);
Zred = zeros(size(Zt,1)-length(ground_elements),size(Zt,1)-length(ground_elements),length(freq));
Yred = Zred;
Tz = zeros(2*size(Zred,1),2*size(Zred,1),length(freq));
Ty = Tz;
Ttotal = Ty;
Ct = Yred;
Gt = Yred;
for a = 1:length(freq)
    Rtemp = NaN(size(R,1));
    Rtemp = real(Z(a))/2.*ones(size(Rtemp,1));
    Rtemp(logical(eye(size(R,1)))) = real(Z(a));
    R(:,:,a) = Rtemp;
    Ltemp = ones(size(L,1)).*imag(Z(a))./w(a);
    Lt(:,:,a) = L+Ltemp;
    
    Zt(:,:,a) = R(:,:,a)+1j*w(a)*Lt(:,:,a);
    red_matrix = eye(size(Zt,1));
    red_matrix(:,ground_elements) =[];
    Zred(:,:,a) = inv(red_matrix'/(Zt(:,:,a))*red_matrix);
    Zred(:,:,a) = (Zred(:,:,a)+transpose(Zred(:,:,a)))/2;
    Tz(:,:,a) = [eye(size(Zred,1)),tsv_length*Zred(:,:,a); zeros(size(Zred,1)), eye(size(Zred,1))];
    
    
    
    
    Ytemp1 = G+1j*w(a)*C;
    Ytemp2 = 1j*w(a)*Cnn;
    Ytemp =  Ytemp1;

    diag_elements = logical(eye(size(Ytemp,1)));
    offdiag_elements = logical((ones(size(Ytemp,1)))-diag_elements);
    Ztemp = inv(Ytemp);
    Ztemp(:) = 1/Ytemp2;
    Ztemp(diag_elements) = Ztemp(diag_elements)+1/Ytemp2;
    Ztemp = Ztemp+inv(Ytemp1);
    Ytemp = inv(Ztemp);

    
    Ytemp_red = red_matrix'*Ytemp*red_matrix;
    Yt(:,:,a) = Ytemp;
    Yred(:,:,a) = Ytemp_red;
    Yred(:,:,a) = (Yred(:,:,a)+transpose(Yred(:,:,a)))/2;

    Ct(:,:,a) = imag(Ytemp_red)./w(a);
    Gt(:,:,a) = real(Ytemp_red);
    Ty(:,:,a) = [eye(size(Zred,1)),zeros(size(Zred,1));tsv_length*Yred(:,:,a) , eye(size(Zred,1))];
    
    
    
    c1 = issymmetric(Zred(:,:,a));
    c2 = issymmetric(Yred(:,:,a));
    
    if (c1 ~= 1 || c2 ~= 1)
        disp('check');
    end
    R1(:,:,a) = real(Zred(:,:,a));
    Lt1(:,:,a) = imag(Zred(:,:,a))./w(a);
    Gt1(:,:,a) = real(Yred(:,:,a));
    Ct1(:,:,a) = imag(Zred(:,:,a))./w(a);
    Ttotal(:,:,a) = Tz(:,:,a)*Ty(:,:,a);
    Zconv(:,:,a) = Ttotal(1:size(Zred,1),size(Zred,1)+1:2*size(Zred,1),a);
    
end


T_TSV = Ttotal;


