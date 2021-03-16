clc
clear
%

% load data
load neuron_data
[nNeuron nSamp] = size(bigdata);
% 
nK    = 7;
Ex    = 0.1*rand(nK,nNeuron);
Bx    = 0.1*rand(nSamp,nK);
%
pix    = exp(-Bx*Ex);pix = pix./(1+pix);xs = bigdata';
%
nPr = nK*(nSamp + nNeuron);
eta = 0.0005;
filen = strcat('TMI_',num2str(nK),'.mat');
%load(filen);
flg = 0;
while flg == 0
    % pi
    pix = exp(-Bx*Ex); 
    pix = pix./(1+pix);pix(isnan(pix)) = 1;
    delt = (pix - xs);
    
    % gradient
    dE = (Bx'*delt);
    dB = ((delt*Ex')); 
    % update
    Ex = Ex + eta*dE;
    Bx = Bx + eta*dB;
    
    % book keeping
    m1 = min(min(pix));m2 = 1-max(max(pix));
    if m1 < 1e-10 || m2 < 1e-10
        flg = 1;
    end
end
save(filen,'Bx','Ex')

%
