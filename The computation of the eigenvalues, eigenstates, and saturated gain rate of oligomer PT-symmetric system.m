clear;
close all;

% set displacement
displacement = 0:5:50;

% components of the source resonator
L_s = 520e-09;
Rs_s = 0.3;
C_s = 216.5e-12;

% components of the neutral resonators
L_n = 520e-09;
Rs_n = 0.3;
C_n = 216.5e-12;

% components of the load resonator
L_l = 520e-09;
Rs_l = 0.3;
C_l = 216.5e-12;
R = 2000;

% calculate resonant frequency and loss rate for efficiency
omega0 = 1/sqrt(L_l*C_l);
gamma_l = 1/(2*R*C_l);

Q_s = omega0*L_s/Rs_s;
Rp_s = Rs_s*(Q_s^2+1);
gamma_s0 = 1/(2*Rp_s*C_s); 

Q_n = omega0*L_n/Rs_n;
Rp_n = Rs_n*(Q_n^2+1);
gamma_n0 = 1/(2*Rp_n*C_n); 
 
Q_l = omega0*L_l/Rs_l;
Rp_l = Rs_l*(Q_l^2+1);
gamma_l0 = 1/(2*Rp_l*C_l);

gamma_lp = gamma_l0+gamma_l;

% set coupling coefficients between neutral resonator and load resonator
CsList = [168.5 178.5 177.5 180.5 168.5 181.5 190.5 216.5 216.5 216.5 216.5 ];

mu1List = [-0.0200971772 -0.0152460508 -0.011129791 -0.0075923041 -0.00494105402 -0.00354020904 -0.00262798069 -0.00202020628 -0.0015885704 -0.00121875547 -0.000975370585 ];
mu2List = [0.2458691880 0.1805948490 0.064527268 -0.00550013859 -0.0200971772 -0.0152460508 -0.011129791 -0.0075923041 -0.0049410540 -0.00354020904 -0.00262798069 ];
mu3List = [-0.0190971772 -0.0055001386 0.0645272680 0.1805948490 0.2458691880 0.1805948490 0.0665272680  -0.0055001386 -0.0200971772 -0.0152460508 -0.0111297910 ];

CsList = CsList.*1e-12;

% set coupling coefficients between neutral resonators
c12 = -0.0610415161;
c23 = -0.0610415161;
c31 = -0.0059353675;

%% set kappa and mu using the flip
mu1 = zeros(1,length(mu1List));
mu2 = zeros(1,length(mu1List));
mu3 = zeros(1,length(mu1List));
kappa12 = zeros(1,length(mu1List));
kappa23 = zeros(1,length(mu1List));
kappa31 = zeros(1,length(mu1List));

for muIdx = 1:length(mu1List)
    u1 = mu1List(muIdx);
    u2 = mu2List(muIdx);
    u3 = mu3List(muIdx);
    
    c12 = -0.0610415161;
    c23 = -0.0610415161;
    c31 = -0.0059353675;
    
    
    if u1 < 0
        u1 = -u1;
        c12 = -c12;
        c31 = -c31;
    end
    if u2 < 0
        u2 = -u2;
        c12 = -c12;
        c23 = -c23;
    end
    if u3 < 0
        u3 = -u3;
        c31 = -c31;
        c23 = -c23;
    end

    mu1(muIdx) = (omega0*u1)/(2*(1-u1^2));
    mu2(muIdx) = (omega0*u2)/(2*(1-u2^2));
    mu3(muIdx) = (omega0*u3)/(2*(1-u3^2));
    
    kappa12(muIdx) = (omega0*c12)/(2*(1-c12^2));
    kappa23(muIdx) = (omega0*c23)/(2*(1-c23^2));
    kappa31(muIdx) = (omega0*c31)/(2*(1-c31^2));
end


%% determine PT-symmetric region and PT-asymmetric region according to displacement
domega = zeros(4, 1);
omegadelta = zeros(4,length(mu1List));
ngsat = zeros(4, length(mu1List));
gsat = zeros(4, length(mu1List));
PTsymmetricIdx = [];
asymmIdx = [];
flag = 0;
for muIdx = 1:length(mu1List)
    kappa = zeros(3,1);
    mu = zeros(3,1);

    kappa(1) = kappa23(muIdx)./omega0;
    kappa(2) = kappa31(muIdx)./omega0;
    kappa(3) = kappa12(muIdx)./omega0;
    
    mu(1) = mu1(muIdx)./omega0;
    mu(2) = mu2(muIdx)./omega0;
    mu(3) = mu3(muIdx)./omega0;
    
    ngamma_lp = gamma_lp/omega0;
    resonantShift = 1;

    % imaginary part of determinant of H
    poly = [1, ...
            (-sum(mu)), ...
            -(sum(kappa.^2)+sum(mu.^2)) + resonantShift * sum(prod(mu)./mu),...
            -2*sum(prod(mu)./mu.*kappa) - 2*prod(kappa) + ...
            (sum(mu)*sum(mu.^2) -sum(mu.^3) + sum(mu.*kappa.^2) - prod(mu)),...
            -2*sum(prod(mu.*kappa)./mu./kappa) + sum(mu.^2.*kappa.^2) + ...
            (-prod(mu)*sum(mu) + 2*prod(mu)*sum(kappa))];
    
    % domega = omega0 - omega       
    domega = sort(roots(poly),'descend');
    omegadelta(:,muIdx) = domega;
    
    % real part of determinant of H
    numerator = domega.^5 +...
              resonantShift * (-sum(mu) .* domega.^4) +  ...
              (-2*sum(mu.^2)-sum(kappa.^2)) .* domega.^3 + ...
              resonantShift * sum((prod(mu)./mu)) * domega.^3 +...
              (-2*prod(kappa) - 4*sum(prod(mu)./mu.*kappa)) .* domega.^2 + ...
              resonantShift * (-prod(mu) + sum(mu.*kappa.^2) + ...
              2*sum(prod(mu)./mu.*(sum(mu)-mu))) * domega.^2 + ...
              (2.*sum(mu.^2.*kappa.^2) - 4.*sum(prod(mu).*prod(kappa)./(mu.*kappa))).* domega + ...
              resonantShift * (prod(mu)*(-2*sum(mu)+4*sum(kappa)) * domega);
    
    denominator = -ngamma_lp*domega.^3 +...
                  resonantShift * ngamma_lp * sum(mu) * domega.^2 + ...
                  ngamma_lp*sum(kappa.^2)*domega +...
                  resonantShift * ngamma_lp * sum(-prod(mu)./mu) * domega +...
                  ngamma_lp*2*prod(kappa) +...
                  resonantShift * ngamma_lp * (sum(-mu.*kappa.^2) + prod(mu));
                      
    % normalized gsat according to branch
    ngsat(:,muIdx) = numerator./denominator/ngamma_lp;

    if min(ngsat(:,muIdx))>= 1
        PTsymmetricIdx = [PTsymmetricIdx, muIdx];
        flag = 1;
    elseif flag == 0
        PTsymmetricIdx = [PTsymmetricIdx, muIdx];
    elseif flag == 1
        asymmIdx = [asymmIdx, muIdx];
    end
end

% simulaiton result
ADS_freq = [11.9 12.58 13.93 12.63 11.95 12.62 13.92 14.87 14.71 14.75 14.77];

ADSa1_mag = [1.717 1.682 1.594 1.684 1.703 1.683 1.602 1.765 1.697 1.772 1.775];
ADSa2_mag = [0.561 0.703 1.812 0.019 0.016 0.026 0.6 0.447 0.912 0.726 0.548];
ADSa3_mag = [2.926 3.001 3.484 0.622 0.586 0.713 2.003 0.082 0.557 0.353 0.238];
ADSa4_mag = [0.561 0.604 0.206 2.949 2.937 2.933 3.394 0.646 0.666 0.627 0.484];
ADSa5_mag = [1.768 1.694 1.598 1.686 1.73 1.694 1.596 0.1 0.33 0.227 0.136];

ADSVec = [ADSa2_mag./ADSa1_mag; ADSa3_mag./ADSa1_mag;...
ADSa4_mag./ADSa1_mag;ADSa5_mag./ADSa1_mag];

% PT-symmetric region
eigenRealFreqList = zeros(length(mu1List),5);
eigenImagFreqList = zeros(length(mu1List),5);
eigenVecList = zeros(length(mu1List),5,5);
optVecList = zeros(length(mu1List),4);
branchList = zeros(length(mu1List),1);

for muIdx = 1:length(mu1List)
    
    L = [[1i*omega0+gamma_lp,   -1i*mu1(muIdx),                           -1i*mu2(muIdx),                           -1i*mu3(muIdx),                                 0];... 
        [-1i*mu1(muIdx),         1i*(omega0-resonantShift*mu1(muIdx)),    -1i*kappa12(muIdx),                       -1i*kappa31(muIdx),                       -1i*mu1(muIdx)];...
        [-1i*mu2(muIdx),        -1i*kappa12(muIdx),                        1i*(omega0-resonantShift*mu2(muIdx)),    -1i*kappa23(muIdx),                       -1i*mu2(muIdx)];...
        [-1i*mu3(muIdx),        -1i*kappa31(muIdx),                       -1i*kappa23(muIdx),                        1i*(omega0-resonantShift*mu3(muIdx)),    -1i*mu3(muIdx)];...
        [0,                     -1i*mu1(muIdx),                           -1i*mu2(muIdx),                           -1i*mu3(muIdx),                            1i*omega0-gamma_lp]];
    
    H = L/(1i);
    [Evector, Evalue] = eig(H,'vector');
    omegaE = Evalue;
    Evalue = Evalue/(2*pi);
    [RealEvalueSorted,sortIdx] = sort(real(Evalue));
    eigenRealFreqList(muIdx,:) = RealEvalueSorted;
    eigenImagFreqList(muIdx,:) = imag(Evalue(sortIdx));
    eigenVecList(muIdx,:,1) = Evector(:,sortIdx(1));
    eigenVecList(muIdx,:,2) = Evector(:,sortIdx(2));
    eigenVecList(muIdx,:,3) = Evector(:,sortIdx(3));
    eigenVecList(muIdx,:,4) = Evector(:,sortIdx(4));
    eigenVecList(muIdx,:,5) = Evector(:,sortIdx(5));
    
    calB1 = eigenVecList(muIdx,2:5,1)/eigenVecList(muIdx,1,1);
    calB2 = eigenVecList(muIdx,2:5,2)/eigenVecList(muIdx,1,2);
    calB3 = eigenVecList(muIdx,2:5,3)/eigenVecList(muIdx,1,3);
    calB4 = eigenVecList(muIdx,2:5,4)/eigenVecList(muIdx,1,4);
    calB5 = eigenVecList(muIdx,2:5,5)/eigenVecList(muIdx,1,5);
    
    tempADS = ADSVec(:,muIdx);
    adsBs = transpose(repmat(tempADS,1,5));
    calBs = [calB1; calB2; calB3; calB4; calB5];
    cost = sum(abs(adsBs-abs(calBs)),2);
    [~,branchIdx] = min(cost);
    branchList(muIdx) = branchIdx;
    optVecList(muIdx,:) = calBs(branchIdx,:);

end

% PT-asymmetric region
brokenBranchList = zeros(length(mu1List),4);
brokenFreqList = zeros(length(mu1List),4);
for muIdx = 1:length(mu1List)
   
    brokenFreqList(muIdx,:) = transpose(sort(1./2./pi*omega0.*(1+omegadelta(:,muIdx))));
    
    [~,lowestgsatIdx] = min(ngsat(:,muIdx));
    tempFreq = brokenFreqList(muIdx,lowestgsatIdx);
    tempOmega = 2*pi*tempFreq;

        A = [[1i*(omega0-tempOmega-mu1(muIdx)),      -1i*kappa12(muIdx),                       -1i*kappa31(muIdx),                       -1i*mu1(muIdx)];...
            [-1i*kappa12(muIdx) ,                      1i*(omega0-tempOmega-mu2(muIdx)),       -1i*kappa23(muIdx),                       -1i*mu2(muIdx)];...
            [-1i*kappa31(muIdx),                      -1i*kappa23(muIdx),                        1i*(omega0-tempOmega-mu3(muIdx)),       -1i*mu3(muIdx)];...
            [-1i*mu1(muIdx),                         -1i*mu2(muIdx),                          -1i*mu3(muIdx),                           1i*(omega0-tempOmega)-gamma_lp]];
        b = [1i*mu1(muIdx);                           1i*mu2(muIdx);                           1i*mu3(muIdx);                           0];
        a = A\b;
        brokenBranchList(muIdx,:) = a;
end

% select EigenFrequency and EigenVector for PT-symmetric and PT-asymmetric region
Result = zeros(length(mu1List),4);
Imag_part = zeros(length(mu1List),5);
FreqList = zeros(length(mu1List),5);

Result(PTsymmetricIdx,:) = abs(optVecList(PTsymmetricIdx,:));
Result(asymmIdx,:) = abs(brokenBranchList(asymmIdx,:));
FreqList(PTsymmetricIdx,:) = eigenRealFreqList(PTsymmetricIdx,:);
FreqList(asymmIdx,1) = brokenFreqList(asymmIdx,1);
FreqList(asymmIdx,2) = brokenFreqList(asymmIdx,2);
FreqList(asymmIdx,3) = brokenFreqList(asymmIdx,3);
FreqList(asymmIdx,4) = brokenFreqList(asymmIdx,3);
FreqList(asymmIdx,5) = brokenFreqList(asymmIdx,4);
Imag_part(PTsymmetricIdx,:) = eigenImagFreqList(PTsymmetricIdx,:);
Imag_part(muIdx,1:4) = imag(brokenFreqList(muIdx,:));

% measurement result
cList_meas = [176.5 172.5 159.5 163.5 161.5 162.5 185 202.5 198.5 198.5 198.5 ];

freq_meas = [11.847 12.415 13.71 12.667 12.02 12.656 13.7 14.677 14.82 14.8 14.82 ];

V1_meas = [4.52 4.52 4.36 4.56 4.6 4.52 4.4 4.6 4.64 4.6 4.6 ];
V2_meas = [1.44 1.76 3.6 0.4 0.4 0.32 0.64 1.52 2.32 2.4 2.56 ];
V3_meas = [7.68 8.08 8.32 2.08 1.6 2.08 3.92 1.28 1.44 1.52 1.52 ];
V4_meas = [1.6 2 2.16 8.72 8.24 8.72 8.96 0.88 1.92 1.92 1.92 ];
V5_meas = [4.48 4.52 4.3 4.52 4.56 4.48 4.36 0.36 0.56 0.48 0.4 ];

cList_meas = cList_meas.*1e-12;

a2_meas = V2_meas./V1_meas;
a3_meas = V3_meas./V1_meas;
a4_meas = V4_meas./V1_meas;
a5_meas = V5_meas./V1_meas;

% calculation efficiency
efficiency_sim = zeros(length(displacement),1);
efficiency_meas = zeros(length(displacement),1);

for muIdx = 1 : length(displacement)
    efficiency_sim(muIdx,1) = (gamma_l*C_l*ADSa5_mag(1,muIdx)^2)/((CsList(muIdx)*gamma_s0*ADSa1_mag(1,muIdx)^2)+(C_s*gamma_n0*ADSa2_mag(1,muIdx)^2)+(C_s*gamma_n0*ADSa3_mag(1,muIdx)^2)+(C_s*gamma_n0*ADSa4_mag(1,muIdx)^2)+(C_s*gamma_lp*ADSa5_mag(1,muIdx)^2));
    efficiency_meas(muIdx,1) = (C_l*gamma_l*V5_meas(1,muIdx)^2)/((cList_meas(muIdx)*gamma_s0*V1_meas(1,muIdx)^2)+(C_s*gamma_n0*V2_meas(1,muIdx)^2)+(C_s*gamma_n0*V3_meas(1,muIdx)^2)+(C_s*gamma_n0*V4_meas(1,muIdx)^2)+(C_s*gamma_lp*V5_meas(1,muIdx)^2));
end

% plot parameter
width = 3.0;
height = 4.0;
lw = 2;
alw = 0.75;
fsz = 10;

% real eigenvalue plot
figure(1);
pos = get(gcf, 'Position');
plot(displacement,FreqList/1e6,'Linewidth',lw); hold on;
plot(displacement,ADS_freq,'--','Color','black','Linewidth',lw); hold on;
plot(displacement,freq_meas,'*','Color','r'); hold on;
xlabel('Displacement [mm]');
ylabel('Frequency [MHz]');
set(gcf, 'Position', [pos(1)-200*width pos(2) width*100, height*100]);
set(gca, 'FontSize', fsz, 'LineWidth', alw);
ax = gca;
ax.FontSize = 11;
xlim([0 displacement(end)]);
ylim([10 22]);
legend('Cal 1.','Cal 2.','Cal 3.','Cal 4.','Cal 5.','Sim.','Mea.');

% imaginary eigenvalue plot
figure(2);
pos = get(gcf, 'Position');
plot(displacement,(Imag_part)/1e6,'Linewidth',lw);
xlabel('Displacement [mm]');
ylabel('Frequency [MHz]');
set(gcf, 'Position', [pos(1)-100*width pos(2) width*100, height*100]);
set(gca, 'FontSize', fsz, 'LineWidth', alw);
legend({'Cal 1.','Cal 2.','Cal 3.','Cal 4.','Cal 5.'},'Location','northwest');
xlim([0 displacement(end)]);
ylim([-0.5 0.5]);

% efficiency plot
figure(3);
pos = get(gcf, 'Position');
plot(displacement,efficiency_sim,'Linewidth',lw,'Color','black'); hold on;
plot(displacement,efficiency_meas,'Linewidth',lw,'Color','red'); hold on;
xlabel('Displacement [mm]');
ylabel('Efficiency');
ax = gca;
ax.FontSize = 11;
set(gcf, 'Position', [pos(1)-100*width pos(2) width*100, height*100]);
legend('Sim.','Mea.');
xlim([0 displacement(end)]);
ylim([0 1]);

% branch gsat plot
ref = ones(1,length(mu1List));
figure(4);
pos = get(gcf, 'Position');
semilogy(displacement,ngsat(1,:),'Linewidth',lw,'Color',[0 0.4470 0.7410]); hold on;
semilogy(displacement,ngsat(2,:),'Linewidth',lw,'Color',[0 0.4470 0.7410]); hold on;
semilogy(displacement,ngsat(3,:),'Linewidth',lw,'Color',[0 0.4470 0.7410]); hold on;
semilogy(displacement,ngsat(4,:),'Linewidth',lw,'Color',[0 0.4470 0.7410]); hold on;
plot(displacement,ref,'Linewidth',lw,'Color',[0 0.4470 0.7410]);
xlabel('Displacement [mm]');
ylabel('g_s/\gamma^{\prime}_l');
text(17,0.6,'g_s/\gamma^{\prime}_l=1');
set(gcf, 'Position', [pos(1)-100*width pos(2) width*100, height*100]);
xlim([15 displacement(end)]);
ax = gca;
ax.FontSize = 12;


width = 12.0;
height = 4.0;
lw = 2;
alw = 0.75;
fsz = 10;

% eigenstate plot
figure(5);
set(gcf, 'Position', [pos(1) pos(2) width*100, height*100]);
subplot(1,4,1);
pos = get(gcf, 'Position');
plot(displacement, Result(:,1),'Linewidth',lw); hold on;
plot(displacement,ADSVec(1,:),'--','Linewidth',lw,'Color','black'); hold on;
plot(displacement,a2_meas,'*','Color','red');
xlabel('Displacement [mm]');
ylabel('V_1/V_s');
set(gca, 'FontSize', fsz, 'LineWidth', alw);
legend('Cal.','Sim.','Mea.');
ax = gca;
ax.FontSize = 11;
xlim([0 displacement(end)]);
ylim([0 3]);

subplot(1,4,2);
pos = get(gcf, 'Position');
plot(displacement, Result(:,2),'Linewidth',lw); hold on;
plot(displacement,ADSVec(2,:),'--','Linewidth',lw,'Color','black'); hold on;
plot(displacement,a3_meas,'*','Color','red');
xlabel('Displacement [mm]');
ylabel('V_2/V_s');
set(gcf, 'Position', [pos(1) pos(2) width*100, height*100]);
set(gca, 'FontSize', fsz, 'LineWidth', alw);
legend('Cal.','Sim.','Mea.');
ax = gca;
ax.FontSize = 11;
xlim([0 displacement(end)]);
ylim([0 3]);

subplot(1,4,3);
pos = get(gcf, 'Position');
plot(displacement, Result(:,3),'Linewidth',lw); hold on;
plot(displacement,ADSVec(3,:),'--','Linewidth',lw,'Color','black'); hold on;
plot(displacement,a4_meas,'*','Color','red');
xlabel('Displacement [mm]');
ylabel('V_3/V_s');
set(gcf, 'Position', [pos(1) pos(2) width*100, height*100]);
set(gca, 'FontSize', fsz, 'LineWidth', alw);
legend('Cal.','Sim.','Mea.');
ax = gca;
ax.FontSize = 11;
xlim([0 displacement(end)]);
ylim([0 3]);


subplot(1,4,4);
pos = get(gcf, 'Position');
plot(displacement, Result(:,4),'Linewidth',lw); hold on;
plot(displacement,ADSVec(4,:),'--','Linewidth',lw,'Color','black'); hold on;
plot(displacement,a5_meas,'*','Color','red');xlabel('Displacement [mm]');
ylabel('V_l/V_s');
set(gcf, 'Position', [pos(1) pos(2) width*100, height*100]);
set(gca, 'FontSize', fsz, 'LineWidth', alw);
legend('Cal.','Sim.','Mea.');
ax = gca;
ax.FontSize = 11;
xlim([0 displacement(end)]);
ylim([0 3]);



