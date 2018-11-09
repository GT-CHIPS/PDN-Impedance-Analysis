function Tgrid = calc_T_grid(himd,cond_si,esi,epoly,cond_cu,tand_poly,t_metal,ws,spacing1,M,N,in_port,out_port,freq,nvia)
%Copyright (c) 2018 Hakki Mert Torun
%Power Delivery Network (PDN) Impedance Analysis for Heteregenous
%Integration. Components in PDN analysis include PCB P/G plane, C4 bump
%array, TSV array, via array, u-bump array.
%This material is based on work supported by DARPA CHIPS project under
%Award N00014-17-1-2950.
%For questions and queries, please contact: htorun3@gatech.edu


% Z-Matrix calculation of Power/Ground grid in interposer based on 
% pi-equivalent unit cells.
% Note that Z-Matrix is converted to T-Matrix.

% The implementation is based on Transmission Matrix method
% Kim, Joong-Ho, and Madhavan Swaminathan. 
%"Modeling of irregular shaped power distribution planes
% using transmission matrix method." 
% IEEE Transactions on Advanced Packaging 24.3 (2001): 334-346.

%Unit cell values are as in:
%Kim, Jaemin, et al. "Chip-package hierarchical power distribution network 
%modeling and analysis based on a segmentation method." 
%IEEE Transactions on advanced packaging 33.3 (2010): 647-659.
%%
mu0 = 4*pi*1e-7;
e0 = 8.854*1e-12;

w = 2*pi*freq;

Rs = 1./(cond_cu*t_metal);
%%
R = 0.25*Rs.*spacing1./ws;
L = 1e-12*(1e6*spacing1.*(0.13.*exp(-spacing1*1e6/45)+0.14*log(spacing1./ws)+0.07));

Ci = epoly*1e-3*((44-28*himd*1e6).*(ws*1e6).^2 + (280*himd*1e6 + 0.8*spacing1*1e6-64).*(ws*1e6)+12*spacing1*1e6-1500*himd+1700);
Cf_nom1 = 4*spacing1.*ws.*1e12.*(log(spacing1*1e6)-log(spacing1*1e6-2*ws*1e6)+exp(-1/3));
Cf_denom1 = ws*1e6*pi + 2*himd*1e6.*(log(spacing1*1e6)-log(spacing1*1e6-2*ws*1e6)+exp(-1/3));
Cf = 1e-15*e0*epoly*1e9*(Cf_nom1./Cf_denom1 + 2*spacing1*1e6/pi*sqrt(himd./(spacing1*1e6-2*1e6*ws)));
C = Cf + 1e-15*Ci;
k = cos(ws*pi/(2*spacing1));
Csi_temp = e0/2*elliprat(k);
Csi_tot = (spacing1-2*ws)*(2*esi*(1+cond_si./(1j*esi*w))*Csi_temp);
Gsi = imag(Csi_tot); Csi = real(Csi_tot); 
Zs = (R+2j*w.*L);
Yp = (w.*(C.*tand_poly+Gsi)+1j*w.*(C+Csi));


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
Z = Cp;
input_ports2keep = round((in_port-sqrt(nvia)/2+0.5):1:(in_port+sqrt(nvia)/2));
output_ports2keep = round((out_port-sqrt(nvia)/2+0.5):1:(out_port+sqrt(nvia)/2));

for a = 1:sqrt(nvia)-1
    input_ports2keep = [input_ports2keep,input_ports2keep(end-sqrt(nvia)+1:end)+M];
    output_ports2keep = [output_ports2keep,output_ports2keep(end-sqrt(nvia)+1:end)+M];
end

for a = 1:length(freq)
    Cp(:,:,a) = full(blktridiag(YMM_inner(:,:,a),-1*YMN(:,:,a),-1*YMN(:,:,a),M));
    Cp(1:N,1:N,a) = YMM(:,:,a);
    Cp(M*N-N+1:M*N,M*N-N+1:M*N,a) = YMM(:,:,a);
    Z(:,:,a) = inv(Cp(:,:,a));

end

Z2 = [Z(in_port,in_port,:) Z(in_port,out_port,:);Z(out_port,in_port,:) Z(out_port,out_port,:)];

Tgrid = z2abcd(Z2);




function rat = elliprat(k)
    kp = sqrt(1-k^2);
    if(k^2<0.5)
        rat = pi/(log(2*(1+sqrt(kp^2))./(1-sqrt(kp^2))));
    else
        rat = 1/pi*(log(2*(1+sqrt(k^2))./(1-sqrt(k^2))));
    end
end
end