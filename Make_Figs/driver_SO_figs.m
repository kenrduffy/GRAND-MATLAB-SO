% Soft output figures plotted output from simulations with both
% ORBGRAND (soft detection), ORBGRAND1 (soft detection). 

clear

n_CODES = 0;

n=32;
k=26;
code.class = 'RLC';

DECODER='GRAND';
n_CODES=n_CODES+1;
filename = ['../RESULTS/SO_' DECODER '_' code.class '_' num2str(n) '_' num2str(k) '_1.mat'];
load(filename,'code');
code.LT = '-d';
code.color = 'k';
codes(n_CODES).code = code; 

DECODER='ORBGRAND';
n_CODES=n_CODES+1;
filename = ['../RESULTS/SO_' DECODER '_' code.class '_' num2str(n) '_' num2str(k) '_1.mat'];
load(filename,'code');
code.LT = '-d';
code.color = 'b';
code.color = [0.1250, 0.6940, 0.9290];
codes(n_CODES).code = code; 

DECODER='ORBGRAND1';
n_CODES=n_CODES+1;
filename = ['../RESULTS/SO_' DECODER '_' code.class '_' num2str(n) '_' num2str(k) '_1.mat'];
load(filename,'code');
code.LT = '-d';
code.color = 'b';
codes(n_CODES).code = code; 

DECODER='SGRAND';
n_CODES=n_CODES+1;
filename = ['../RESULTS/SO_' DECODER '_' code.class '_' num2str(n) '_' num2str(k) '_1.mat'];
load(filename,'code');
code.LT = '-d';
code.color = 'r';
codes(n_CODES).code = code;

make_fig(codes,1)

%%
function make_fig(codes,fig_no)

    FONT=16;
    MS=8;
    LW=2;

    % Number of conditional probability groups
    n_groups = 50;

    % For each code
    for ii=1:length(codes)
        % SO
        SO = codes(ii).code.DEC_OUTPUT(:,1);
        % Decoding correctness
        dec_correct = codes(ii).code.DEC_OUTPUT(:,2);

        x=zeros(1,n_groups);
        f=zeros(1,n_groups);
        n_elements = size(codes(ii).code.DEC_OUTPUT,1);
        N = floor(n_elements/n_groups);
        for jj=1:n_groups
            these = find((jj-1)*1/n_groups<SO & SO<=jj*1/n_groups);
            x(jj) = mean(SO(these));
            f(jj) = mean(dec_correct(these));  
        end

        codes(ii).code.x = x;
        codes(ii).code.f = f;
    end

    figure(fig_no)
    clf
    hold on
    for ii=1:length(codes)
            code_info = [codes(ii).code.DECODER ', ' codes(ii).code.class ' [' num2str(codes(ii).code.n) ',' num2str(codes(ii).code.k) '], R=' num2str(codes(ii).code.k/codes(ii).code.n,'%.2f')];
            plot(codes(ii).code.x,codes(ii).code.f,'-d','displayname',code_info,'color',codes(ii).code.color,'LineWidth',LW, 'MarkerSize',MS)
    end
    plot([0 1],[0 1],'--','displayname','x=y','LineWidth',LW,'MarkerSize',MS)
    hold off
    ylabel('Empirical conditional probability correct')
    xlabel('Predicted probability correct')
    grid 'on'
    legend('show','Location','SouthEast');
    set(gca,'FontSize',FONT)
    title('Forecaster Calibration')

 
end
