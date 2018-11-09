function T_bump = calc_T_bump(cond_cu,rbump,pitch,bump_length,freq)
%%
%Copyright (c) 2018 Hakki Mert Torun
%Power Delivery Network (PDN) Impedance Analysis for Heteregenous
%Integration. Components in PDN analysis include PCB P/G plane, C4 bump
%array, TSV array, via array, u-bump array.
%This material is based on work supported by DARPA CHIPS project under
%Award N00014-17-1-2950.
%For questions and queries, please contact: htorun3@gatech.edu

% C4-bumps and u-bumps are approximated as cynlindrical structures.
% Multiconductor RLGC matrices of these are calculated as same as in
% "calc_T_TSV" script but without depletion capacitance. See the comments
% in there for explanations.



%%
mu0 = 4*pi*1e-7;
e0 = 8.854*1e-12;
esi = 11.9;
epoly = 3.9;



%%
rvia = rbump;
tox = 0e-6;
tdep = 0e-6;

%%

columnNumber = 10;
rowNumber = 10;

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
%%
Lii = mu0/2/pi*2*log(dn0/(rvia+tox+tdep));
L = dn0'*dn0;
L = mu0/2/pi*log(L/(rvia+tox+tdep));
L = L-mu0/2/pi*log(dij);
L(1:rowNumber*columnNumber:end) = Lii;

%%

%%

C = mu0*e0*inv(L);

G = 0;

Lii = mu0/2/pi*2*log(dn0/(rvia));
L = dn0'*dn0;
L = mu0/2/pi*log(L/(rvia));
L = L-mu0/2/pi*log(dij);
L(1:rowNumber*columnNumber:end) = Lii;

%%
Cnn = (log((rvia+tox)/rvia)/(2*pi*epoly*e0)+log((rvia+tox+tdep)/(rvia+tox))/(2*pi*esi*e0))^-1;

w = 2*pi*freq;
Z = sqrt(1j*w*mu0/cond_cu)/(2*pi*(rvia)).*besseli(0,rvia*sqrt(1j*w*mu0*cond_cu))./besseli(1,rvia*sqrt(1j*w*mu0*cond_cu));
Lt = zeros(size(L,1),size(L,2),length(freq));
R = Lt;

Zt = Lt; 
Yt = Lt; 

%%
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
    Tz(:,:,a) = [eye(size(Zred,1)),bump_length*Zred(:,:,a); zeros(size(Zred,1)), eye(size(Zred,1))];
    
    
    
    
    Ytemp1 = G+1j*w(a)*C;
    Ytemp =  Ytemp1;


    Ytemp_red = red_matrix'*Ytemp*red_matrix;
    Yt(:,:,a) = Ytemp;
    Yred(:,:,a) = Ytemp_red;
    Yred(:,:,a) = (Yred(:,:,a)+transpose(Yred(:,:,a)))/2;

    Ct(:,:,a) = imag(Ytemp_red)./w(a);
    Gt(:,:,a) = real(Ytemp_red);
    Ty(:,:,a) = [eye(size(Zred,1)),zeros(size(Zred,1));bump_length*Yred(:,:,a) , eye(size(Zred,1))];
    
    
    
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
    
    
end

T_bump = Ttotal;
