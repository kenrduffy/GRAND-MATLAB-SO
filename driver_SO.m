% Ken R. Duffy, 2018-2025.

% Guessing Random Additive Noise Decoding (GRAND)
% All code is subject to license:
% GRAND Codebase Non-Commercial Academic Research Use License 021722.pdf

% This performs decoding and collates soft output (SO) forecasts for decoding
% correctness.

% GRAND algorithms efficiently decode any code of any structure when the 
% number of redunant bits, n-k, is moderate. 

% The algorithms are highly parallizable, which results in low latency when
% full utilized, as described in in-silicon implementations for soft detection 
% ORBGRAND:
% A. Riaz, A. Yasar, F. Ercan, W. An, J. Ngo, K. Galligan, M. Médard, 
% K. R. Duffy & R. T. Yazicigil. "A sub-0.8 pJ/bit, 16.3 Gbps/mm2 
% universal soft-detection decoder using ORBGRAND in 40 nm CMOS." 
% Proceedings of IEEE ISSCC, 2023.

% The inefficient MATLAB implementations here, however, do not exploit the 
% parallelizability and so obtaining the full performance of codes with 
% n-k>20 may be time-consuming.

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
% 70, 4528–4542, 2022 and implemented here using landlside and the version 
% introduced in 
% K. Galligan, M. Médard, K. R. Duffy, "Block turbo decoding with ORBGRAND"
% Conference on Information Science and Systems, 2023.

% Soft outputfrom soft input decoders follows from: 
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

% Decoders are: 
% 'GRAND' (hard detection decoding with soft detection APP)
% 'ORBGRAND' (soft detection); 
% 'ORBGRAND1' (soft detection);
% 'SGRAND' (soft detection);
DECODER ='GRAND';

% Sim range in Eb/N0
ebn0=-2:0.5:8;

%GRAND parameters
max_query = inf; % Can reduce to an abandonment value.

% Create this many decodings at each Eb/N0 
N_decodings = 10^3;

% Modulation schemes available using MATLAB's toolbox:
modlist = {'pi/2-BPSK','BPSK','QPSK','16QAM','64QAM','256QAM'}; 
bpsList = [1 1 2 4 6 8];
% Pick the modulation
modulation = 'BPSK';
% Determine the number of bits per symbol in the modulation
nmodbits = bpsList(strcmpi(modlist,modulation));

% Pick the code
% RLC, eRLC (even RLC), PAC, CAPOLAR, BCH, eBCH or CRC. 
code.class = 'RLC';

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
    disp(['With n-k=' num2str(n-k) ' expect slow decoding in noisy channels with this inefficient MATLAB implementation.'])
end

% Check to see if it's an even code, which can be exploited in ORBGRAND to
% reduce complexity
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

% For each decoding, record the SO, whether the decoding is correct, the
% number of queries until a codeword is found, and the ebn0
DEC_OUTPUT = [];
total_decoded = 0;

for ii=1:length(snr_db)
    % To measure BLER in real time
    zz=0;
    
    % Noise variance
    sigma2 = 1/(10^(0.1*snr_db(ii))); 
    % Using MATLAB's channel function
    awgnchan = comm.AWGNChannel('NoiseMethod','Variance','Variance',sigma2);

    % Decode N_decodings at each Eb/N0
    while num_decoded(ii)<N_decodings

        num_decoded(ii)=num_decoded(ii)+1;
        total_decoded=total_decoded+1;
                
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
        correct = isequal(y_decoded, c);
     
        DEC_OUTPUT = [DEC_OUTPUT; APP correct n_guess ebn0(ii)];

        % Report to the terminal every time 10^3 decodings performed
        if mod(total_decoded,10^3)==0
            disp(['n=' num2str(n) ', k=' num2str(k) ', R=' num2str(k/n,'%.2f') ', ' num2str(ebn0(ii)) ' dB, Dec=' num2str(num_decoded(ii)) ', Total=' num2str(total_decoded)]); 
        end
        
    end
end


% Save the results.

if isinf(max_query)
    filename = ['../RESULTS/SO_' DECODER '_' code.class '_' num2str(n) '_' num2str(k) '_' num2str(nmodbits) '.mat'];
else
    filename = ['../RESULTS/SO_' DECODER '_' code.class '_' num2str(n) '_' num2str(k) '_' num2str(nmodbits) '_' num2str(log2(max_query)) '.mat'];
end


% Store the data in a way where results from future simulations with the 
% same code can be appended.

if exist(filename, 'file') == 2
    load(filename,'code')
    code.DEC_OUTPUT = [code.DEC_OUTPUT; DEC_OUTPUT];
else
    code.n = n;
    code.k = k;
    code.modulation = modulation;
    code.nmodbits = nmodbits;
    code.G = G;
    code.DEC_OUTPUT = DEC_OUTPUT;
    code.max_query = max_query;
    code.DECODER=DECODER;
end

snr_db = unique(code.DEC_OUTPUT(:,4));
code.snr_db = snr_db;

n_snr = length(snr_db);
BLER = zeros(1,n_snr); % Block error rate
BS = zeros(1,n_snr); % Brier score

for ii=1:n_snr
    these = find(code.DEC_OUTPUT(:,4)==snr_db(ii));
    BLER(ii) = sum(code.DEC_OUTPUT(these,2)==0)/length(these);
    BS(ii)= sum((code.DEC_OUTPUT(these,2)-code.DEC_OUTPUT(these,1)).^2)/length(these);
end
code.BLER = BLER;
code.BS = BS;

save(filename,'code')

disp(['Code=' code.class ', ' num2str(sum(num_decoded)) ' packets decoded.'])




