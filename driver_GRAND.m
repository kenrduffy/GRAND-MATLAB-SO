% Ken R. Duffy, 2018-2025.

% Guessing Random Additive Noise Decoding (GRAND)
% All code is subject to license:
% GRAND Codebase Non-Commercial Academic Research Use License 021722.pdf

% This runs the main simulation of GRAND/ORBGRAND/ORBGRAND1/SGRAND.

% GRAND algorithms efficiently decode any code of any structure when the 
% number of redunant bits, n-k, is moderate. In the presence of Soft Input,
% GRAND algorithms can readily produce accurate blockwise Soft Output.

% The algorithms are highly parallizable, which results in low latency when
% full utilized, as described in in-silicon implementations for hard 
% detection GRAND, e.g.:
% A. Riaz, V. Bansal, A. Solomon, W. An, Q. Liu, K. Galligan, K. R. Duffy, 
% M. Médard & R. T. Yazicigil. "Multi-code multi-rate universal maximum 
% likelihood decoder using GRAND", Proceedings of IEEE ESSCIRC, 2021.
% and soft detection ORBGRAND:
% A. Riaz, A. Yasar, F. Ercan, W. An, J. Ngo, K. Galligan, M. Médard, 
% K. R. Duffy & R. T. Yazicigil. "A sub-0.8 pJ/bit, 16.3 Gbps/mm2 
% universal soft-detection decoder using ORBGRAND in 40 nm CMOS." 
% Proceedings of IEEE ISSCC, 2023.

% The inefficient MATLAB implementations here are designed for illustration.
% They do not exploit the parallelizability and so obtaining the full 
% performance of codes with n-k>20 may be time-consuming. That is
% particularly the case with SGRAND.

% GRAND
% K. R. Duffy, J. Li, and M. Medard, "Capacity-achieving guessing random 
% additive noise decoding," IEEE Trans. Inf. Theory, vol. 65, no. 7, pp. 
% 4023–4040, 2019.

% Ordered Reliability Bits GRAND (ORBGRAND)
% Basic ORBGRAND as introduced in
% K. R. Duffy, “Ordered reliability bits guessing random additive noise 
% decoding," in IEEE ICASSP, 2021, pp. 8268–8272 
% and implemented using the Landslide algorithm introduced in 
% K. R. Duffy, W. An, and M. Medard, "Ordered reliability bits guessing 
% random additive noise decoding,” IEEE Transactions on Signal Processing, 
% 70, 4528–4542, 2022.

% 1-line ORBGRAND as introduced in
% K. R. Duffy, W. An, and M. Medard, "Ordered reliability bits guessing 
% random additive noise decoding," IEEE Transactions on Signal Processing, 
% 70, 4528–4542, 2022 and implemented here using the landlside algorithm 
% and the integer intercept determination described in 
% K. Galligan, M. Médard, K. R. Duffy, "Block turbo decoding with ORBGRAND"
% Conference on Information Science and Systems, 2023.

% Soft Input GRAND (SGRAND)
% A. Solomon, K. R. Duffy and M. Medard, "Soft maximum likelihood decoding 
% using GRAND", IEEE ICC, 2020.

% Soft output from soft input decoders follows from: 
% K. Galligan, P. Yuan, M. Médard & K. R. Duffy. Upgrade error detection 
% to prediction with GRAND". Proceedings of Globecom, 2023.
% P. Yuan, M. Médard, K. Galligan & K. R. Duffy, "Soft-output (SO) GRAND 
% and long, low rate codes to outperform 5 LDPC codes", IEEE Trans. 
% Wireless Commun., 24(4), 3386-3399, 2025.

% Exploiting codebook structure to discard query patters in GRAND
% algorithms was first proposed in 
% M. Rowshan and J. Yuan, "Constrained Error Pattern Generation for GRAND", 
% IEEE International Symposium on Information Theory, 2022. 
% Availing of the fact that the implementation of GRAND, ORBGRAND,
% ORBGRAND1 means that patterns with certain Hamming weights 
% can be skipped without creation, here we exploit Even Codes, e.g. 
% Lin & Costello, "Error control coding: fundamentals and applications", 
% 2004. This results in fewer queries but identical performance. Further
% query reduction is possible with additional methods.

% When codebook structure is used to inform query order, the Soft Output
% computation is adapted. Here, that is applied as in
% J. Feng, K. R. Duffy and M. Medard, "Leveraging Code Structure to Improve 
% Soft Output for GRAND, GCD, OSD, and SCL", arXiv:2503.16677.

clear;

% Decoder is: 
% 'GRAND' (hard detection); 
% 'ORBGRAND' (soft detection); 
% 'ORBGRAND1' (soft detection);
% 'SGRAND' (soft detection); 

DECODER='ORBGRAND1';
% Sim range in Eb/N0
ebn0=-2:0.5:8;

%GRAND parameters
max_query = inf; % Can reduce to an abandonment value.

% Go to this many errors for each Eb/N0 value
err_thresh = 50;

% Modulation schemes available using MATLAB's toolbox:
modlist = {'pi/2-BPSK','BPSK','QPSK','16QAM','64QAM','256QAM'}; 
bpsList = [1 1 2 4 6 8];
% Pick the modulation
modulation = 'BPSK';
% Determine the number of bits per symbol in the modulation
nmodbits = bpsList(strcmpi(modlist,modulation));

% Pick the code
% RLC, eRLC (even RLC), POLAR, CAPOLAR, BCH, eBCH or CRC. 
code.class = 'CRC';

% Code dimensions
n=32;
k=26;

% Create the generation and parity check matrices and record the generator
% polynomial if a CRC
[G,H,hex_poly]=make_code(code.class,n,k);
code.G=G;
code.H=H;
% Record the polynomial
code.poly = hex_poly;

if (n-k)>20
    disp(['With n-k=' num2str(n-k) ' expect slow decoding in noisy channels with this MATLAB implementation.'])
end

% Check to see if it's an even code, which can be exploited in GRAND, 
% ORBGRAND and ORBGRAND 1, but not SGRAND, to reduce query numbers by only
% generating queries whose parity is the same as the demodulated sequence
even_code = 0;
if isequal(mod(sum(G,2),2),zeros(k,1))
    even_code = 1;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Convert Eb/N0 to SNR
snr_db = ebn0+10*log10(k/n)+10*log10(nmodbits);

% Pad modulation with meaningless bits if necessary.
mod_pad=zeros(1,mod(n,nmodbits));

% For each SNR record the following
num_decoded = zeros(1,length(snr_db)); % Number of packets decoded
num_demod_errs = num_decoded; % Number of demodulations that are in error
num_demod_bit_errs = num_decoded; % Number of erroneous bits
num_errs = num_decoded; % Number of erroneous decodings
num_bit_errs = num_decoded; % Number of erroneous bits
num_aband = num_decoded; % Number of abandoned decodings
num_queries = num_decoded; % Total number of code-book queries
tot_app = num_decoded; % To enable average APP at each SNR

for ii=1:length(snr_db)
    tic
    % To measure BLER in real time
    zz=0;
    
    % Noise variance
    sigma2 = 1/(10^(0.1*snr_db(ii))); 
    % Using MATLAB's channel function
    awgnchan = comm.AWGNChannel('NoiseMethod','Variance','Variance',sigma2);

    % Keep decoding packets until you see err_thresh erroneous decodings
    % while num_errs(ii)<err_thresh
    while min(num_decoded(ii)-num_errs(ii),num_errs(ii))<err_thresh

        num_decoded(ii)=num_decoded(ii)+1;
                
        % Encode
        % Uniform random information word
        u = binornd(ones(1,k),0.5);
        % Codeword
        c = mod(u*G,2);
	    % Modulate
	    modOut = nrSymbolModulate([c mod_pad]',modulation);
        % Add White Gaussian noise
        rSig = awgnchan(modOut);
        % Soft demodulate
        y_soft = nrSymbolDemodulate(rSig,modulation,sigma2);
    	y_soft = y_soft(1:n)';
        % Hard demodulate
        y_demod = (y_soft<0);

        % Count the times the demodulation is in error
        if ~isequal(y_demod, c)
            num_demod_errs(ii) = num_demod_errs(ii)+1;
            % Count the demodulated bit errors
            num_demod_bit_errs(ii) = num_demod_bit_errs(ii)+sum(abs(y_demod-c));
        end
        
        % Decode with GRAND (hard detection) but calculates SO if optional
        % y_soft SI given 
        if isequal(DECODER,'GRAND')
            [y_decoded,~,n_guess,abandoned,APP] = bin_GRAND(H,max_query,y_demod,even_code,y_soft);    
        % Decode with basic ORBGRAND (soft detection)
        elseif isequal(DECODER,'ORBGRAND')
            [y_decoded,~,n_guess,abandoned,APP] = bin_ORBGRAND(H,max_query,y_soft,even_code); 
        % Decode with 1-line ORBGRAND (soft detection)
        elseif isequal(DECODER,'ORBGRAND1')
            [y_decoded,~,n_guess,abandoned,APP] = bin_ORBGRAND1(H,max_query,y_soft,even_code); 
        % Decode with SGRAND (soft detection)
        elseif isequal(DECODER,'SGRAND')
            [y_decoded,~,n_guess,abandoned,APP] = bin_SGRAND(H,max_query,y_soft); 
        else
            disp('DECODER must be GRAND or ORBGRAND or ORBGRAND1 or SGRAND')
            return;
        end

        % Total number of queries made at this SNR
        num_queries(ii) = num_queries(ii) + n_guess;
        % Enable averag APP evaluation
        tot_app(ii) = tot_app(ii) + APP;

        % If there is an error in the decoding
        if ~isequal(y_decoded, c)
            % Increment the number of errors observed at this SNR
            num_errs(ii)=num_errs(ii)+1;
            % Count the bit errors.
            if ~abandoned
                % Extract information bits from non-systematic codes
                if isequal(code.class,'CAPOLAR')
                    u_decoded = gflineq(G',single(y_decoded'),2);
                    u_decoded = u_decoded';
                else
                    u_decoded = y_decoded(1:k);
                end
                num_bit_errs(ii) = num_bit_errs(ii) + sum(abs(u_decoded-u));
            else
                % If we've abandoned, then use the demodulated sequence 
                % to determine how many bit errors there were.
                num_bit_errs(ii) = num_bit_errs(ii) + sum(abs(y_demod-c));
                num_aband(ii) = num_aband(ii)+1;
            end
            % Report to the terminal every time a new error is observed
            if mod(num_errs(ii),10)==0
                disp(['n=' num2str(n) ', k=' num2str(k) ', R=' num2str(k/n,'%.2f') ', ' num2str(ebn0(ii)) ' dB, Dec=' num2str(num_decoded(ii)) ', Err=' num2str(num_errs(ii)) ', BLER=' num2str(num_errs(ii)/num_decoded(ii)) ', BER=' num2str(num_bit_errs(ii)/(num_decoded(ii)*k)) ', EG=' num2str(num_queries(ii)/num_decoded(ii),'%.0f')]); 
            end
        end
    end
    toc
end

% Save the results.

if isinf(max_query)
    filename = ['RESULTS/' DECODER '_' code.class '_' num2str(n) '_' num2str(k) '_' num2str(nmodbits) '.mat'];
else
    filename = ['RESULTS/' DECODER '_' code.class '_' num2str(n) '_' num2str(k) '_' num2str(nmodbits) '_' num2str(log2(max_query)) '.mat'];   
end

% Store the data in a way where results from future simulations with the 
% same code can be appended.

if exist(filename, 'file') == 2
    load(filename,'code')
    for ii=1:length(snr_db)
        these = find(code.snr==snr_db(ii));
        if isempty(these)
            code.snr(end+1) = snr_db(ii);
            code.num_decoded(end+1)=num_decoded(ii);
            code.num_demod_errs(end+1)=num_demod_errs(ii);
            code.num_demod_bit_errs(end+1)=num_demod_bit_errs(ii);
            code.num_errs(end+1)=num_errs(ii);
            code.num_bit_errs(end+1)=num_bit_errs(ii);
            code.num_aband(end+1)=num_aband(ii);
            code.num_queries(end+1)=num_queries(ii);
            code.tot_app(end+1)=tot_app(ii);
        else
            code.num_decoded(these)=code.num_decoded(these)+num_decoded(ii);
            code.num_demod_errs(these)=code.num_demod_errs(these)+num_demod_errs(ii);
            code.num_demod_bit_errs(these)=code.num_demod_bit_errs(these)+num_demod_bit_errs(ii);
            code.num_errs(these)=code.num_errs(these)+num_errs(ii);
            code.num_bit_errs(these)=code.num_bit_errs(these)+num_bit_errs(ii);
            code.num_aband(these)=code.num_aband(these)+num_aband(ii);
            code.num_queries(these)=code.num_queries(these)+num_queries(ii);
            code.tot_app(these)=code.tot_app(these)+tot_app(ii);
        end
    end     
else
    code.n = n;
    code.k = k;
    code.modulation = modulation;
    code.nmodbits = nmodbits;
    code.G = G;
    code.snr = snr_db;
    code.num_decoded=num_decoded;
    code.num_demod_errs=num_demod_errs;
    code.num_demod_bit_errs=num_demod_bit_errs;
    code.num_errs=num_errs;
    code.num_bit_errs=num_bit_errs;
    code.num_aband=num_aband;
    code.num_queries=num_queries;
    code.tot_app=tot_app;
end

% Sort according to increasing snr;
[~,ord]=sort(code.snr);
code.snr = code.snr(ord);
code.num_decoded = code.num_decoded(ord);
code.num_demod_errs = code.num_demod_errs(ord);
code.num_demod_bit_errs = code.num_demod_bit_errs(ord);
code.num_errs = code.num_errs(ord);
code.num_bit_errs = code.num_bit_errs(ord);
code.num_aband = code.num_aband(ord);
code.num_queries = code.num_queries(ord);
code.tot_app = code.tot_app(ord);

% Summary statistics   
code.ebn0 = code.snr-10*log10(k/n)-10*log10(nmodbits); %Eb/N0
code.R = code.k/code.n; % Code rate
code.BLERdemod = code.num_demod_errs./code.num_decoded; % Demod BLER
code.BERdemod = code.num_demod_bit_errs./(code.num_decoded*code.n); % Demod BER
code.BLER = code.num_errs./code.num_decoded; % Decoded BLER
code.BER = code.num_bit_errs./(code.num_decoded*code.k); % Decoded BER
code.BLAB= code.num_aband./code.num_decoded; % Decoded Abandoment Rate
code.BLERnoab = (code.num_errs-code.num_aband)./code.num_decoded; % Decoded BLER without abandonment
code.BLERab = (code.num_errs-code.num_aband)./(code.num_decoded-code.num_aband); % Nonabandoned BLER
code.EG = code.num_queries./code.num_decoded; % Average number of code-book queries
code.APP = code.tot_app./code.num_decoded; % Average a posteriori probability of correctness
code.max_query = max_query;

save(filename,'code')

disp(['Decoder=' DECODER ', code=' code.class ', ' num2str(sum(num_decoded)) ' packets decoded.'])




