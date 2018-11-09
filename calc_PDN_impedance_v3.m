%%
%Copyright (c) 2018 Hakki Mert Torun
%Power Delivery Network (PDN) Impedance Analysis for Heteregenous
%Integration. Components in PDN analysis include PCB P/G plane, C4 bump
%array, TSV array, via array, u-bump array.
%This material is based on work supported by DARPA CHIPS project under
%Award N00014-17-1-2950.
%For questions and queries, please contact: htorun3@gatech.edu


clear all
close all
clc

%% Material Properties
cond_si = 10; %S/m
esi = 11.9;
epoly = 3.9;
cond_cu = 5.96e7;
tand_poly = 0.002;
tand_PCB = 0.018;
%% Geometrical Properties
hsub_PCB = 2.54e-5; %d1
epcb = 4;
t_metal_PCB = 36e-6;
unit_cell_size_PCB = 5e-3;
%Number of unit cells of PCB P/G plane(M -> width, N -> Height)
M_PCB = 10;
N_PCB = 10;
%Input and Output Ports in MxN grid
in_port_PCB = 25;
out_port_PCB = 85;
%% P/G Grid Parameters
himd = 1e-6;
t_metal = 1e-6;
grid_width = 40e-6;
grid_spacing = 100e-6;
M_grid = 50;
N_grid = 50;
in_port_grid = M_grid+N_grid/2;
out_port_grid = in_port_grid+4*M_grid;
%% TSV array Parameters
tsv_length = 200e-6;
tsv_radius = 5e-6;
tsv_pitch = 100e-6;
tsv_oxide_thickness = 0.1e-6;
%% C4 bump array Parameters
C4_radius = 50e-6;
C4_pitch = 100e-6;
C4_height = 2*C4_radius;

%% u-bump array Parameters
ubump_radius = 12.125e-6;
ubump_pitch = 45e-6;
ubump_height = 2*ubump_radius;

%% via array parameters
via_radius = 0.5e-6;
via_pitch = 1.5e-6;
via_height = 2*himd+2*0.4e-6; %This is for 5 layers
DC_via = via_height./(cond_cu*pi*via_radius.^2);
%% VRM Inductance can be included here
VRM_inductance = 0;
%% Calculate T-Matrix of Each Component.
%FREQUENCY RANGE THE PDN IMPEDANCE IS CALCULATED
freq = logspace(-3,1,200)*1e9; 
Z_VRM = 1j*2*pi*freq*VRM_inductance;
for a = 1:length(freq)
    T_VRM(:,:,a) = [1, Z_VRM(a); 0 1];
end
Z_PCB_Plane = calc_Z_plane(hsub_PCB,epcb,cond_cu,tand_PCB,t_metal_PCB,unit_cell_size_PCB,M_PCB,N_PCB,in_port_PCB,out_port_PCB,freq);
T_TSV       = calc_T_TSV(cond_si,esi,epoly,cond_cu,tsv_length,tsv_radius,tsv_oxide_thickness,tsv_pitch,freq);
T_C4        = calc_T_bump(cond_cu,C4_radius,C4_pitch,C4_height,freq);
T_ubump     = calc_T_bump(cond_cu,ubump_radius,ubump_pitch,ubump_height,freq);
T_via       = calc_T_bump(cond_cu,via_radius,via_pitch,via_height,freq);

tic
%This takes around 5-6 minutes to calculate due to large M & N 
%(to be improved)
T_grid      = calc_Z_grid(himd,cond_si,esi,epoly,cond_cu,tand_poly,t_metal,grid_width,grid_spacing,M_grid,N_grid,in_port_grid,out_port_grid,freq,25);
% Load existing T-grid if available
% load('T_grid_himd_1um_gridwidth_40um_spacing_100um_M_50_N_50_in_75_out_25.mat');
toc
%% Calculate Contribution of each component to PDN impedance
Y_TSV = abcd2y(T_TSV);
Y_ubump = abcd2y(T_ubump);
Y_C4 = abcd2y(T_ubump);
Y_via = abcd2y(T_via);
Y_grid = abcd2y(T_grid);
Y_PCB_Plane = z2y(Z_PCB_Plane);

Y_TSV_SA = Y_TSV(25,25,:);
Z_TSV_SA = squeeze(1./Y_TSV_SA);

Y_ubump_SA = Y_ubump(25,25,:);
Z_ubump_SA = squeeze(1./Y_ubump_SA);

Y_via_SA = Y_via(25,25,:);
Z_via_SA = squeeze(1./Y_via_SA);

Y_C4_SA = Y_C4(25,25,:);
Z_C4_SA = squeeze(1./Y_C4_SA);

Y_grid_SA = Y_grid(1,1,:);
Z_grid_SA = squeeze(1./Y_grid_SA);

Y_PCB_SA = Y_PCB_Plane(1,1,:);
Z_PCB_SA = squeeze(1./Y_PCB_SA);
%% Cascade here
T_PCB_Plane = z2abcd(Z_PCB_Plane);
% T_TSV = z2abcd(Z_TSV);
% T_C4 = z2abcd(Z_C4);
% T_ubump = z2abcd(Z_ubump);
% T_via = z2abcd(Z_via);
% T_grid = z2abcd(Z_grid);

T_PDN = NaN(2,2,length(freq));
Z_PDN = NaN(1,length(freq));
Z_Ref = 50;

%Connector Matrix for parallel excitation of bumps
connector_matrix = zeros(size(T_via,1),2);
connector_matrix(1:size(T_via,1)/2,1) = 1;
connector_matrix(size(T_via,1)/2+1:end,2) = 1;

%Set "only_ivr" to false for VRM on PCB + IVR on interposer
%Set "only_ivr" to true for only IVR on interposer
only_ivr = false;

%Matrix reduction for parallel excitation
%Based on: Chapter 4 of 
%Young, Brian. Digital signal integrity: modeling and simulation with 
%interconnects and packages. Prentice Hall PTR, 2000.

for a = 1:length(freq)
    
    if (only_ivr)
        T_temp1 = T_ubump(:,:,a)*T_via(:,:,a);
    else
        T_temp1 = T_C4(:,:,a)*T_TSV(:,:,a);
    end
    T_temp2 = T_via(:,:,a)*T_ubump(:,:,a);
    Y_temp1 = abcd2y(T_temp1);
    Y_temp2 = abcd2y(T_temp2);
    
    Y_temp11 = connector_matrix'*Y_temp1*connector_matrix;
    Y_temp22 = connector_matrix'*Y_temp2*connector_matrix;
    
    T_TX = y2abcd(Y_temp11);
    T_RX = y2abcd(Y_temp22);
    
    if (only_ivr)
        T_PCB(:,:,a) = eye(size(T_TX,1));
    else
        T_PCB(:,:,a) = T_VRM(:,:,a)*T_PCB_Plane(:,:,a)*T_TX;
    end
    T_INTERP(:,:,a) = T_grid(:,:,a)*T_RX;
   
    T_PDN(:,:,a) = T_PCB(:,:,a)*T_grid(:,:,a)*T_RX;
    
end

%Short Port 1 to get PDN Impedance, Z_PDN
Y_PDN_FULL = abcd2y(T_PDN);
Y_PDN = Y_PDN_FULL(2,2,:);
Z_PDN = squeeze(1./Y_PDN);

%% Plot PDN Impedance and Breakdown of Component impedances
figure
h = gca;
loglog(freq,abs(Z_PDN),'LineWidth',7)
hold on
loglog(freq,abs(Z_PCB_SA),'LineWidth',3)
loglog(freq,abs(Z_C4_SA),'LineWidth',3)
loglog(freq,abs(Z_TSV_SA),'LineWidth',3)
loglog(freq,abs(Z_via_SA),'LineWidth',3)
loglog(freq,abs(Z_ubump_SA),'LineWidth',3);
loglog(freq,abs(Z_grid_SA),'LineWidth',3);
legend('Z_{PDN}','Z_{PCB}','Z_{C4}','Z_{TSV}','Z_{via}','Z_{ubump}','Z_{grid}')

xlabel('Frequency')
ylabel('PDN Impedance [\Omega]')
h.XGrid = 'on'; h.XMinorGrid = 'on';
h.YGrid = 'on'; h.YMinorGrid = 'on';
h.FontWeight = 'bold'; h.FontSize = 16;
h.XLim = [0.001,10]*1e9;
