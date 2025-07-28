
clear

code.class = 'eBCH';
n=32;
k=26;

n_GRAND = 0;
n_GRAND = n_GRAND+1;
DECODER='GRAND';

filename = strcat(['../RESULTS/' DECODER '_' code.class '_' num2str(n) '_',num2str(k) '_1.mat']);
load(filename,'code');
code.decoder = DECODER;
code.LT = '-d';
code.color='k';
codes(n_GRAND).code = code;

n_GRAND = n_GRAND+1;
DECODER='ORBGRAND';
filename = strcat(['../RESULTS/' DECODER '_' code.class '_' num2str(n) '_',num2str(k) '_1.mat']);
load(filename,'code');
code.decoder = DECODER;
code.LT = '-d';
code.color='b';
codes(n_GRAND).code = code;

n_GRAND = n_GRAND+1;
DECODER='ORBGRAND1';
filename = strcat(['../RESULTS/' DECODER '_' code.class '_' num2str(n) '_',num2str(k) '_1.mat']);
load(filename,'code');
code.decoder = DECODER;
code.LT = '-d';
code.color='g';
codes(n_GRAND).code = code;


n_GRAND = n_GRAND+1;
DECODER='SGRAND';
filename = strcat(['../RESULTS/' DECODER '_' code.class '_' num2str(n) '_',num2str(k) '_1.mat']);
load(filename,'code');
code.decoder = DECODER;
code.LT = '-d';
code.color='r';
codes(n_GRAND).code = code;


make_fig(codes,1)

%%
clear

n=128;
k=113;

n_GRAND = 0;


n_GRAND = n_GRAND+1;
DECODER='ORBGRAND1';
code.class = 'eRLC';
filename = strcat(['../RESULTS/' DECODER '_' code.class '_' num2str(n) '_',num2str(k) '_1.mat']);
load(filename,'code');
code.decoder = DECODER;
code.LT = '-d';
code.color='k';
codes(n_GRAND).code = code;

n_GRAND = n_GRAND+1;
DECODER='ORBGRAND1';
code.class = 'CAPOLAR';
filename = strcat(['../RESULTS/' DECODER '_' code.class '_' num2str(n) '_',num2str(k) '_1.mat']);
load(filename,'code');
code.decoder = DECODER;
code.LT = '-d';
code.color='b';
codes(n_GRAND).code = code;

n_GRAND = n_GRAND+1;
DECODER='ORBGRAND1';
code.class = 'CRC';
filename = strcat(['../RESULTS/' DECODER '_' code.class '_' num2str(n) '_',num2str(k) '_1.mat']);
load(filename,'code');
code.decoder = DECODER;
code.LT = '-d';
code.color='g';
codes(n_GRAND).code = code;

n_GRAND = n_GRAND+1;
DECODER='ORBGRAND1';
code.class = 'eBCH';
filename = strcat(['../RESULTS/' DECODER '_' code.class '_' num2str(n) '_',num2str(k) '_1.mat']);
load(filename,'code');
code.decoder = DECODER;
code.LT = '-d';
code.color='r';
codes(n_GRAND).code = code;


make_fig(codes,2)



%%

function make_fig(codes,fig_no)

    FONT=16;
    MS=8;
    LW=2;
     
    num_codes = length(codes);

    
    figure(fig_no)
    clf
    hold on
    for ii=1:num_codes
        plot(codes(ii).code.ebn0,codes(ii).code.BLER,codes(ii).code.LT,'displayname',[codes(ii).code.decoder ', ' codes(ii).code.class ' [' num2str(codes(ii).code.n) ',' num2str(codes(ii).code.k) '], R=' num2str(codes(ii).code.k/codes(ii).code.n,'%.2f')],'color',codes(ii).code.color,'LineWidth',LW, 'MarkerSize',MS)
        plot(codes(ii).code.ebn0,codes(ii).code.BER,'--d','HandleVisibility','off','color',codes(ii).code.color,'LineWidth',LW, 'MarkerSize',MS)
    end
    hold off
    legend('show','Location','SouthWest');
    ylabel('BLER (solid) / BER (dashed) ')
    xlabel('Eb/N0')
    grid 'on'
    xl = xlim;
    xlim([xl(1)-1 xl(2)+1])
    set(gca, 'YScale', 'log')
    set(gca,'FontSize',FONT)
  
    
end
