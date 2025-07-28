% Guessing Random Additive Noise Decoding (GRAND)
% All code is subject to license
% GRAND Codebase Non-Commercial Academic Research Use License 021722.pdf

% GRAND (hard detection, BSC)
% As introduced in
% K. R. Duffy, J. Li, and M. Medard, "Capacity-achieving guessing random 
% additive noise decoding," IEEE Trans. Inf. Theory, vol. 65, no. 7, pp. 
% 4023–4040, 2019.

% The decoding does not use the soft information. If soft information is
% provided, it is solely used to create the a posteriori estimate that the 
% decoding is correct.

% Inputs:
%   n               - code length
%   H               - Parity check matrix or CRC function
%   max_query       - Maximum number of code-book queries to abandonment
%   y_demod         - Hard decision bits
%   even_code       - 1 if code is even
%   y_soft (opt)    - Optioncal channel soft information causes soft output
%                     to be calculated

%
% Outputs:
%   y_decoded       - Decoded codeword
%   putative_noise  - noise
%   n_guesses       - Number of guesses performed
%   abandoned       - 1 if abandoned, 0 if found a codeword
%   APP             - soft output, likelihood decoding is correct

function [y_decoded,putative_noise,n_guesses,abandoned,APP] = bin_GRAND(H,max_query,y_demod,even_code,y_soft)

    % Hard demodulation is input to the decoder
    demod_parity = mod(sum(y_demod),2);

    % If there's soft input, calculate SO
    SI= exist('y_soft', 'var');

    [k,n]=size(H);
    k=n-k;
    n_guesses = 0;
    abandoned = 0;

    if SI==1
        % Bit reliability
        abs_LLR = abs(y_soft);
        % Prior belief that the demodulated bit is in error
        p_demod = exp(-abs_LLR)./(1+exp(-abs_LLR));
        % Probability normalizer (all zeros guess)
        pG1 = prod(1-p_demod);
        % Sum of probabilities of guesses
        sPG = 0;
        
        % If the code is even
        if even_code 
            p_even = 1/2*(1+prod(1-2*p_demod));
            if demod_parity==0
                pG1 = pG1/p_even;
            else
                pG1 = pG1/(1-p_even);
            end
        end
    end

    [err_loc_vec, err_vec, ~] = gen_next_err(n); %Generate initial vectors

    % If the code is even, then check the parity of the demod to decide
    % whether to query even or odd
    if even_code
        y_demod_parity = mod(sum(y_demod),2);
        if y_demod_parity==1
            [err_loc_vec, err_vec, ~] = gen_next_err(n, err_loc_vec,1); %Generate next error vector
        end
        HW_inc = 2; % Hamming weight increment
    else
        HW_inc = 1;
    end

    while n_guesses < max_query
        [decoded,putative_noise,ng]=bin_syn_check(H,y_demod,err_vec);
        n_guesses = n_guesses+ng;

        % If there's soft input, calculate soft output while decoding
        if SI==1
            % APP accounting
            if n_guesses==1
                if even_code==0 || y_demod_parity ==0
                    pG=pG1;
                else
                    pG=0;
                end
            else
                pG = pG1*prod(p_demod(err_loc_vec)./(1-p_demod(err_loc_vec))); % likelihood of query
            end
            sPG = sPG+pG;
        end

        if (decoded == 1)
            y_decoded = mod(y_demod-putative_noise,2);

            if SI==1
                if even_code == 0
                    APP = pG/(pG+(1-sPG)*(2^k-1)/(2^n-n_guesses));
                else
                    APP = pG/(pG+(1-sPG)*(2^k-1)/(2^(n-1)-n_guesses));
                end
            else
                APP = NaN;
            end

            return;
        end
        [err_loc_vec, err_vec, ~] = gen_next_err(n, err_loc_vec,HW_inc); %Generate next error vector
    end

    % If abandoned
    y_decoded = -1*ones(size(y_demod));
    abandoned = 1;


end

function [decoded,noise_found,n_guesses]=bin_syn_check(H,y_demod,err_vec)

    decoded = 0;
    noise_found = NaN(1,size(err_vec,2));
    % How many guesses have we made
    n_guesses=0;
    while (decoded == 0 && n_guesses<size(err_vec,1))
        % How much work have we done
        n_guesses=n_guesses+1;
        t = mod(y_demod-err_vec(n_guesses,:),2);
        if (mod(H*t',2) == zeros(size(H,1),1))
            decoded = 1;
            noise_found = err_vec(n_guesses,:);
        end   
    end

end

%
%Description: This function generates the new error location vector,
%given a previous error location vector and mask vector of possible error locations
%
%Inputs:
%   n           - Code length
%   err_loc_vec - Previous error location vector
%
%Outputs:
%   err_loc_vec - New error location vector. [] if a new error location vector cannot be generated
%   err_vec     - Binary error vector that corresponds to err_loc_vec. Zero vector if no error could be generated
%

function [err_loc_vec, err_vec, weight] = gen_next_err(n, err_loc_vec,inc)


    mask_vec = 1:n;
    err_vec = zeros(1,n);

    if ~exist('err_loc_vec', 'var')
        %Initialize zero error location vector
        err_loc_vec = [];
        weight = 0;
        return;
    end


    %Get next error location vector
    err_loc_vec = increase_error(mask_vec, err_loc_vec,inc);
    if isequal(err_loc_vec, []) %Could not generate error
        weight=inf;
        return;
    end

    %Generate error
    weight = length(err_loc_vec);
    err_vec(err_loc_vec) = 1;


    end

    function err_loc_vec = increase_error(mask_vec, err_loc_vec,inc)
    %
    %Description: This function generates the next error location vector given the previous one
    %

    max_possible_err_loc = length(mask_vec);
    %Try to remain in the same weight
    success = false;
    for ii=length(err_loc_vec):-1:1
        new_err_loc_vec = err_loc_vec;
    %     success = true;
        if new_err_loc_vec(ii)==max(mask_vec) %Cannot increase this error
            continue;
        end
        success = true;
        ind_in_mask = find(mask_vec==new_err_loc_vec(ii)); %Find the index in mask_vec such that mask_vec(index)==err_loc_vec(index)
        new_err_loc_vec(ii) = mask_vec(ind_in_mask+1);
        for jj=ii+1:length(err_loc_vec)
            new_ind = ind_in_mask+jj-ii+1;
            if new_ind<=max_possible_err_loc
                new_err_loc_vec(jj) = mask_vec(new_ind);
            else
                success = false;
                break;
            end
        end
        if success==true
            err_loc_vec = new_err_loc_vec;
            break;
        end
    end

    %If need to go to the next weight
    if success==false
        if length(err_loc_vec)==max_possible_err_loc
            err_loc_vec = [];
            return;
        end
        for ii=1:length(err_loc_vec)+inc % inc = 1 or 2
            err_loc_vec(ii) = mask_vec(ii);
        end
    end

end

