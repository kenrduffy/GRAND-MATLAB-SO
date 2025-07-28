% Guessing Random Additive Noise Decoding (GRAND)
% All code is subject to license
% GRAND Codebase Non-Commercial Academic Research Use License 021722.pdf

% SGRAND (soft detection)
% As introduced in
% A. Solomon, K. R. Duffy and M. Medard, "Soft maximum likelihood 
% decoding using GRAND", Proceedings IEEE ICC, 2020.

% Inputs:
%   H               - Parity check matrix
%   max_query    - Number of guesses to perform before abandoning
%   y_soft           - Log Likelihood Ratio log(p(y|0)/p(y|1))
%   full_save_flag  - If true and saving guessing order, will continue guessing regardless of whether a codeword is found. Default: true
%
% Outputs:
%   y_decoded       - Decoded ML codeword
%   putative_noise  - ML noise
%   n_guesses       - Number of guesses performed
%   abandoned       - 1 if abandoned, 0 if found a codeword
%   APP             - soft output, likelihood decoding is correct

function [y_decoded,putative_noise,n_guesses,abandoned,APP] = bin_SGRAND(H,max_query, y_soft)

    max_query = min(1e6,max_query); % Need non-inf for the heap to function.
 
    [k,n]=size(H);
    k=n-k;
    y=(y_soft<0);

    %Default parameters
    n_guesses = 0;
    heap_size = max_query;

    %Abandon defualt values
    % Throw an error if we abandon
    y_decoded = -1*ones(size(y));
    abandoned =1;
    % Bit reliability
    abs_LLR = abs(y_soft);
    % Prior belief that the demodulated bit is in error
    p_demod = exp(-abs_LLR)./(1+exp(-abs_LLR));
    % Probability normalizer (all zeros guess)
    pG1 = prod(1-p_demod);
    % Sum of probabilities of guesses
    sPG = 0;

    
    %Calculate p0_vec, p1_vec - the log pdf assuming the bit is 0', 1' respectively
    lr = exp(y_soft);
    p1_vec = log(1./(1+lr));
    p0_vec = log(1./(1+1./lr));
    [~, ind_order] = sort(min(p0_vec,p1_vec),'descend');

    %Create first error vector - no errors
    err_struct.err_vec = false(1,n);
    err_vec = err_struct.err_vec+0; %Convert to double
    err_struct.score = sum(max(p0_vec,p1_vec));
    %Create Max Heap and add the most likely error to it
    err_heap = MaxHeap(heap_size); 
    err_heap.InsertKey(err_struct);

    
    while n_guesses < max_query
        [err_vec, err_heap, ~] = gen_next_err_soft(err_heap, p0_vec, p1_vec, ind_order); %SGRAND
        [decoded,putative_noise,ng]=bin_syn_check(H,y,err_vec); %Syndrome check
        n_guesses = n_guesses+ng;

        % APP accounting
        if n_guesses==1
            pG=pG1;
        else
            ind=find(err_vec==1);
            pG = pG1*prod(p_demod(ind)./(1-p_demod(ind))); % likelihood of query
        end
        sPG = sPG+pG;

        % If codeword found
        if (decoded == 1)
            y_decoded = mod(y-putative_noise,2);
            abandoned = 0;
            APP = pG/(pG+(1-sPG)*(2^k-1)/(2^n-n_guesses));
            return;
        end
    end
end


function [y] = bpsk_modem(x,mod_demod)
%function [y] = bpsk_modem(x,mod_demod)
%
%Description: This is a modem (modulater demodulater) for BPSK modulation
%               The used convention is 1' -> -1, 0' -> 1
%
%Inputs:
%   x           - Input vector (binary for modulation, real for demodulation)
%   mod_demod   - 'mod' for modulation, 'demod' for demodulation
%
%Outputs:
%   y           - Modulated real vector for modulation, demodulated binary vector for demodulation
%

if strcmp(mod_demod, 'mod')
    y = -(2*x-1); %BPSK modulation: 1' -> -1, 0' -> 1
elseif strcmp(mod_demod, 'demod')
    y = double(x<0); %Demodulation
end

end

% Assumes code is binary & determine code-book membership by syndrome.
% Takes parity check matrix, received signal (vector), and matrix of
% putative noises.
function [decoded,noise_found,n_guesses]=bin_syn_check(H,y,noise)

    decoded = 0;
    noise_found = NaN(1,size(noise,2));
    % How many guesses have we made
    n_guesses=0;
    while (decoded == 0 && n_guesses<size(noise,1))
        % How much work have we done
        n_guesses=n_guesses+1;
        t = mod(y-noise(n_guesses,:),2);
        if (mod(H*t',2) == zeros(size(H,1),1))
            decoded = 1;
            noise_found = noise(n_guesses,:);
        end
    end

end

function [next_err, err_heap, next_err_score] = gen_next_err_soft(err_heap, p0_vec, p1_vec, ind_order)
%function [next_err, err_heap, next_err_score] = gen_next_err_soft(err_heap, p0_vec, p1_vec, ind_order)
%
%Description: This function generates the next error vector to be used, as well as the heap to be used in the next call of this function
%
%Inputs:
%   err_heap        - Errors max heap. If L<inf, then this is a struct containing score_vec and err_struct_vec.
%   p0_vec          - Log probability vector. p0_vec(i) = log(Pr(i-th bit is 0'))
%   p1_vec          - Log probability vector. p1_vec(i) = log(Pr(i-th bit is 1'))
%   ind_order       - Error index order, from most likely error to least likely error (in bits)
%
%Outputs:
%   next_err        - Next error to be checked
%   err_heap        - Errors max heap
%   next_err_score  - Next error score
%

n = length(ind_order);

%Get most likely error
next_err_struct = err_heap.ExtractMax;
err_vec_sorted = next_err_struct.err_vec(ind_order);

last_err_ind = find(err_vec_sorted==true,1, 'last');
if isempty(last_err_ind) %If the error is the 0 error vector, start at the first index
    last_err_ind=0;
end

if last_err_ind < n
    %Add new error with an additional 1 at the next available index
    new_err_struct = next_err_struct;
    new_err_ind = ind_order(last_err_ind+1);
    %Update score
    pmax_new = max(p0_vec(new_err_ind), p1_vec(new_err_ind));
    pmin_new = min(p0_vec(new_err_ind), p1_vec(new_err_ind));
    new_err_struct.score = new_err_struct.score - pmax_new + pmin_new; %update score
    %Flip bit
    new_err_struct.err_vec(new_err_ind) = 1;
    %Update heap
    err_heap.InsertKey(new_err_struct);
%     disp(['Adding: ', num2str(find(new_err_struct.err_vec==1))]); %DEBUG ONLY
    
    if last_err_ind>0 %Bit can be moved
        %Add an error with the same Hamming weight, move the last one
        %Update score
        last_err_ind = ind_order(last_err_ind);
        pmax_old = max(p0_vec(last_err_ind), p1_vec(last_err_ind));
        pmin_old = min(p0_vec(last_err_ind), p1_vec(last_err_ind));
        new_err_struct.score = new_err_struct.score + pmax_old - pmin_old; %update score
        %Unflip current bit, next bit already flipped
        new_err_struct.err_vec(last_err_ind) = 0;
        %Update heap
        err_heap.InsertKey(new_err_struct);
    end
end

next_err = next_err_struct.err_vec+0; %+0 to convert to double
next_err_score = next_err_struct.score;
end

function err_heap = delete_sorted_vec_item(err_heap, ind, replace_flag)
%This function deletes the ind-th item from the sorted errors vector (err_heap)
if replace_flag==true
    %Create dummy element
    err_heap.score_vec = [-inf, err_heap.score_vec];
    err_heap.err_struct_vec = [err_heap.err_struct_vec(1),err_heap.err_struct_vec];
    ind = ind+1;
end
    %Delete item
    err_heap.score_vec(ind) = [];
    err_heap.err_struct_vec(ind) = [];
end


