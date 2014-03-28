% close all
clear all
clc

load eta.dat

load pmf1PDF.dat
load pmf1CV.dat
load pmf1CSDRH.dat
load pmf1CSDRI.dat

load pmf2PDF.dat
load pmf2CV.dat
load pmf2CSDRH.dat
load pmf2CSDRI.dat

load betaPDF.dat
load betaCV.dat
load betaCSDRH.dat
load betaCSDRI.dat

load deltaPDF.dat
load unfrmCV.dat
load unfrmCSDR.dat

load girCSDRH.dat
load amcBetaCSDRH.dat
load amcPmfCSDRH.dat
load linCV.dat

% load PMFIIDERIVS.dat

eta_p = [0.27 0.27];
eta_st = [0.351 0.351];

neta = size(eta,1);

scrsz = get(0,'ScreenSize');
figsize = [2 scrsz(4)/5 scrsz(3)/1.5 scrsz(4)/2];
fntSize = 16;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

h = figure('Position',figsize);
semilogy(eta,pmf2PDF,'k-',...
         eta,pmf1PDF,'r-',...  
         eta,betaPDF,'b-',... 
    'LineWidth',1)
yl = ylim;
xl = xlim;
hold on 
plot(eta_p,yl,'m--','LineWidth',1)
hold on 
plot(eta_st,yl,'m--','LineWidth',1)
hold off
xlabel('$$\eta$$','Interpreter','LaTex','FontSize',fntSize)
ylabel('$$\tilde P(\eta)$$','Interpreter','LaTex','FontSize',fntSize)
leg = legend('Trinary PMF{-}PDF',...
             'Binary PMF{-}PDF',...
             '\beta{-}PDF',... 
             'Location','North');
set(leg,'FontNAme','Times','FontSize',fntSize-6)
set(gca,'FontNAme','Times','FontSize',fntSize)
ylim(yl)
xlim(xl)
text(eta_st(1)+0.01,(yl(1)+yl(2))/2,'$\eta_{st}$','Interpreter','LaTex','FontSize',fntSize-2); 
text(eta_p(1)-0.03,(yl(1)+yl(2))/2,'$\eta_{p}$','Interpreter','LaTex','FontSize',fntSize-2);
tightfig(h);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

h = figure('Position',figsize);
plot(eta,pmf2PDF,'k-',...
     eta,pmf1PDF,'r-',...  
     eta,betaPDF,'b-',... 
    'LineWidth',1)
yl = ylim;
xl = xlim;
hold on 
plot(eta_p,yl,'m--','LineWidth',1)
hold on 
plot(eta_st,yl,'m--','LineWidth',1)
hold off
xlabel('$$\eta$$','Interpreter','LaTex','FontSize',fntSize)
ylabel('$$\tilde P(\eta)$$','Interpreter','LaTex','FontSize',fntSize)
leg = legend('Trinary PMF{-}PDF',...
             'Binary PMF{-}PDF',...
             '\beta{-}PDF',... 
             'Location','North');
set(leg,'FontNAme','Times','FontSize',fntSize-6)
set(gca,'FontNAme','Times','FontSize',fntSize)
ylim(yl)
xlim(xl)
text(eta_st(1)+0.01,(yl(1)+yl(2))/2,'$\eta_{st}$','Interpreter','LaTex','FontSize',fntSize-2); 
text(eta_p(1)-0.03,(yl(1)+yl(2))/2,'$\eta_{p}$','Interpreter','LaTex','FontSize',fntSize-2);
tightfig(h);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

h = figure('Position',figsize);
plot(eta(2:neta-1),pmf2CV(2:neta-1,1),'k-',...
     eta(2:neta-1),pmf1CV(2:neta-1,1),'r-',...
     eta(2:neta-1),betaCV(2:neta-1,1),'b-',... 
     eta(2:neta-1),linCV(2:neta-1,1),'g-',... 
     eta(1:neta),unfrmCV(1:neta,1),'c-',...
     'LineWidth',1)
yl = ylim;
xl = xlim;
hold on 
plot(eta_p,yl,'m--','LineWidth',1)
hold on 
plot(eta_st,yl,'m--','LineWidth',1)
hold off
xlabel('$$\eta$$','Interpreter','LaTex','FontSize',fntSize)
ylabel('$$\langle u\vert\eta\rangle \; [\rm{m}/\rm{s}]$$','Interpreter','LaTex','FontSize',fntSize)
leg = legend('PDF-gradient, Trinary PMF{-}PDF',...
             'PDF-gradient, Binary PMF{-}PDF',...
             'PDF-gradient, \beta{-}PDF',... 
             'Linear',...
             'Mean',...
             'Location','North');
set(leg,'FontNAme','Times','FontSize',fntSize-6)
set(gca,'FontNAme','Times','FontSize',fntSize)
ylim(yl)
xlim(xl)
text(eta_st(1)+0.01,(yl(1)+yl(2))/2,'$\eta_{st}$','Interpreter','LaTex','FontSize',fntSize-2); 
text(eta_p(1)-0.03,(yl(1)+yl(2))/2,'$\eta_{p}$','Interpreter','LaTex','FontSize',fntSize-2);
tightfig(h);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

h = figure('Position',figsize);
plot(eta(2:neta-1),pmf2CV(2:neta-1,2),'k-',...
     eta(2:neta-1),pmf1CV(2:neta-1,2),'r-',...
     eta(2:neta-1),betaCV(2:neta-1,2),'b-',...
     eta(2:neta-1),linCV(2:neta-1,2),'g-',... 
     eta(1:neta),unfrmCV(1:neta,2),'c-',...
     'LineWidth',1)
yl = ylim;
xl = xlim;
hold on 
plot(eta_p,yl,'m--','LineWidth',1)
hold on 
plot(eta_st,yl,'m--','LineWidth',1)
hold off
xlabel('$$\eta$$','Interpreter','LaTex','FontSize',fntSize)
ylabel('$$\langle v\vert\eta\rangle \; [\rm{m}/\rm{s}]$$','Interpreter','LaTex','FontSize',fntSize)
leg = legend('PDF-gradient, Trinary PMF{-}PDF',...
             'PDF-gradient, Binary PMF{-}PDF',...
             'PDF-gradient, \beta{-}PDF',... 
             'Linear',...
             'Mean',...
             'Location','North');
set(leg,'FontNAme','Times','FontSize',fntSize-6)
set(gca,'FontNAme','Times','FontSize',fntSize)
ylim(yl)
xlim(xl)
text(eta_st(1)+0.01,(yl(1)+yl(2))/2,'$\eta_{st}$','Interpreter','LaTex','FontSize',fntSize-2); 
text(eta_p(1)-0.03,(yl(1)+yl(2))/2,'$\eta_{p}$','Interpreter','LaTex','FontSize',fntSize-2);
tightfig(h);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

h = figure('Position',figsize);
plot(eta,pmf2CSDRH,'k--',...
     eta,pmf2CSDRI,'k-',... 
     eta,pmf1CSDRH,'r--',...
     eta,pmf1CSDRI,'r-',...  
     eta,betaCSDRH,'b--',...
     eta,betaCSDRI,'b-',... 
     eta,unfrmCSDR,'c-',...
    'LineWidth',1)
yl = ylim;
xl = xlim;
hold on 
plot(eta_p,yl,'m--','LineWidth',1)
hold on 
plot(eta_st,yl,'m--','LineWidth',1)
hold off
xlabel('$$\eta$$','Interpreter','LaTex','FontSize',fntSize)
ylabel('$$\langle \chi\vert\eta\rangle \; [1/\rm{s}]$$','Interpreter','LaTex','FontSize',fntSize)
leg = legend('Trinary PMF - Homogeneous',...
             'Trinary PMF - Inhomogeneous',...
             'Binary PMF - Homogeneous',...
             'Binary PMF - Inhomogeneous',...
             'Beta - Homogeneous',...
             'Beta - Inhomogeneous',...
             'Mean',...
             'Location','NorthWest');
set(leg,'FontNAme','Times','FontSize',fntSize-6)
set(gca,'FontNAme','Times','FontSize',fntSize)
ylim(yl)
xlim(xl)
text(eta_st(1)+0.01,(yl(1)+yl(2))/2,'$\eta_{st}$','Interpreter','LaTex','FontSize',fntSize-2); 
text(eta_p(1)-0.03,(yl(1)+yl(2))/2,'$\eta_{p}$','Interpreter','LaTex','FontSize',fntSize-2);
tightfig(h);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

h = figure('Position',figsize);
plot(eta,pmf1CSDRH,'k-',...
     eta(1:2:neta),amcPmfCSDRH(1:2:neta),'ko',...
     eta,betaCSDRH,'r-',...
     eta(1:2:neta),girCSDRH(1:2:neta),'ro',...
     'LineWidth',1)
yl = ylim;
xl = xlim;
hold on 
plot(eta_p,yl,'m--','LineWidth',1)
hold on 
plot(eta_st,yl,'m--','LineWidth',1)
hold off
xlabel('$$\eta$$','Interpreter','LaTex','FontSize',fntSize)
ylabel('$$\langle \chi\vert\eta\rangle \; [1/\rm{s}]$$','Interpreter','LaTex','FontSize',fntSize)
% leg = legend('Mortensen (Homogeneous), Binary PMF{-}PDF',...
%              'AMC, Binary PMF{-}PDF',...
%              'Mortensen (Homogeneous), \beta{-}PDF',...
%              'Girimaji, \beta{-}PDF',...
%              'Location','NorthWest');
leg = legend(sprintf('Mortensen (Homogeneous),\nBinary PMF{-}PDF'),...
             'AMC, Binary PMF{-}PDF',...
             sprintf('Mortensen (Homogeneous),\n\\beta{-}PDF'),...
             'Girimaji, \beta{-}PDF',...
             'Location','NorthWest');
set(leg,'FontNAme','Times','FontSize',fntSize-6)
set(gca,'FontNAme','Times','FontSize',fntSize)
ylim(yl)
xlim(xl)
text(eta_st(1)+0.01,(yl(1)+yl(2))/2,'$\eta_{st}$','Interpreter','LaTex','FontSize',fntSize-2); 
text(eta_p(1)-0.03,(yl(1)+yl(2))/2,'$\eta_{p}$','Interpreter','LaTex','FontSize',fntSize-2);
tightfig(h);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
unfrmCSDR(1)
trapz(eta,betaCSDRH.*betaPDF)
trapz(eta,betaCSDRI.*betaPDF)
trapz(eta,pmf1CSDRH.*pmf1PDF)
trapz(eta,pmf1CSDRI.*pmf1PDF)
trapz(eta,pmf2CSDRH.*pmf2PDF)
trapz(eta,pmf2CSDRI.*pmf2PDF)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% h = figure('Position',figsize);
% plot(eta(2:neta-1),PMFIIDERIVS(2:neta-1,1),'k-',...
%      eta(2:neta-1),PMFIIDERIVS(2:neta-1,2),'r-',...
%      eta(2:neta-1),PMFIIDERIVS(2:neta-1,3),'b-',...
%      eta(2:neta-1),pmf2PDF(2:neta-1),'g-',...
%      'LineWidth',1)
% hold on 
% y = ylim;
% ylim(y)
% plot(eta_p,y,'c--','LineWidth',1)
% hold off
% xlabel('$$\eta$$','Interpreter','LaTex','FontSize',fntSize)
% leg = legend('$\partial^2 II/\partial\tilde\xi\tilde\xi$',...
%              '$\partial^2 II/\partial\widetilde{\xi''^2}\widetilde{\xi''^2}$',...
%              '$\partial^2 II/\partial\tilde\xi\partial\widetilde{\xi''^2}$',...
%              '$\tilde P(\eta)$',...
%              'Location','North');         
% set(leg,'FontNAme','Times','FontSize',fntSize-6,'Interpreter','Latex')
% set(gca,'FontNAme','Times','FontSize',fntSize)
% tightfig(h);

