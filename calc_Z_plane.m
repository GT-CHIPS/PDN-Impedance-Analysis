function Z_plane = calc_Z_plane(hsub,epoly,cond_cu,tand_poly,t_metal,unit_cell_size,M,N,in_port,out_port,freq)
%Copyright (c) 2018 Hakki Mert Torun
%Power Delivery Network (PDN) Impedance Analysis for Heteregenous
%Integration. Components in PDN analysis include PCB P/G plane, C4 bump
%array, TSV array, via array, u-bump array.
%This material is based on work supported by DARPA CHIPS project under
%Award N00014-17-1-2950.
%For questions and queries, please contact: htorun3@gatech.edu


% Z-Matrix calculation of Power/Ground plane pairs based on pi-equivalent
% unit cells.
% The implementation is based on Transmission Matrix method
% Kim, Joong-Ho, and Madhavan Swaminathan. 
%"Modeling of irregular shaped power distribution planes
% using transmission matrix method." 
% IEEE Transactions on Advanced Packaging 24.3 (2001): 334-346.
%%
hsi = 0; %d2
mu0 = 4*pi*1e-7;
e0 = 8.854*1e-12;


%%
d = hsub+hsi; %distance between P/G planes
freq = freq';
w = 2*pi*freq;


h = unit_cell_size;

%%
Rdc = 2./(cond_cu*t_metal);
Rac = 2*(1+1j)*sqrt(pi*freq*mu0/cond_cu);
L = mu0*d;
C = epoly*e0*h^2/d;
G = w*C*tand_poly;
%%


Zs = Rdc+Rac+1j*w*L;
Yp = 1j*w*C+G;

%%


YMM = NaN(N,N,length(freq));
YMN = NaN(N,N,length(freq));
ytemp = Yp/2+2./Zs;
ytemp2 = Yp/4+1./Zs;
ztemp = 1./Zs;
nOnes = ones(N,1);
for a = 1:length(freq)
    YMM(:,:,a) = diag(ytemp(a).*nOnes,0) + diag(-1/2*ztemp(a).*nOnes(1:N-1),-1) + diag(-1/2*ztemp(a).*nOnes(1:N-1),1);  
    YMM(1,1,a) = ytemp2(a);
    YMM(end,end,a) = ytemp2(a);
    YMN(:,:,a) = diag(ztemp(a).*nOnes,0);
    YMN(1,1,a) = 0.5*ztemp(a);
    YMN(N,N,a) = 0.5*ztemp(a);
end
YMM_inner = YMM*2;
Cp = NaN(M*N,M*N,length(freq));
Tp = NaN(2*M*N,2*M*N,length(freq));
for a = 1:length(freq)
    Cp(:,:,a) = full(blktridiag(YMM_inner(:,:,a),-1*YMN(:,:,a),-1*YMN(:,:,a),M));
    Cp(1:N,1:N,a) = YMM(:,:,a);
    Cp(M*N-N+1:M*N,M*N-N+1:M*N,a) = YMM(:,:,a);
    Tp(:,:,a) = [eye(M*N),zeros(M*N,M*N); Cp(:,:,a),eye(M*N)];
end
% Z = abcd2z(Tp);

Z = y2z(Cp);
Z_plane = [Z(in_port,in_port,:) Z(in_port,out_port,:); Z(out_port,in_port,:), Z(out_port,out_port,:)];
end

