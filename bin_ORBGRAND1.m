% Guessing Random Additive Noise Decoding (GRAND)
% All code is subject to license
% GRAND Codebase Non-Commercial Academic Research Use License 021722.pdf

% 1-line ORGBRAND (soft detection) 
% As introduced in
% K. R. Duffy, “Ordered reliability bits guessing random additive noise 
% decoding," in IEEE ICASSP, 2021, pp. 8268–8272 
% and implemented using the Landslide algorithm introduced in 
% K. R. Duffy, W. An, and M. Medard, "Ordered reliability bits guessing 
% random additive noise decoding,” IEEE Trans. Signal Process., vol. 70, 
% pp. 4528-4542, 2022.
% Soft output from 
% K. Galligan, P. Yuan, M. Médard & K. R. Duffy. Upgrade error detection 
% to prediction with GRAND". Proceedings of Globecom, 2023.
% P. Yuan, M. Médard, K. Galligan & K. R. Duffy, "Soft-output (SO) GRAND 
% and long, low rate codes to outperform 5 LDPC codes", IEEE Trans. 
% Wireless Commun., 24(4), 3386-3399, 2025.

% 1-line ORBGRAND differs from Basic ORBGRAND by considering a (dynamically
% determined) intercept in its statistical model.

% Inputs:
%   n               - code length
%   H               - Parity check matrix or CRC function
%   max_query       - Maximum number of code-book queries to abandonment
%   y_soft          - Channel soft information
%   even_code       - 1 if code is even
%
% Outputs:
%   y_decoded       - Decoded ML codeword
%   err_vec         - Putative noise
%   n_guesses       - Number of guesses performed
%   abandoned       - 1 if abandoned, 0 if found a codeword
%   APP             - soft output, likelihood decoding is correct

function [y_decoded,err_vec,n_guesses,abandoned,APP] = bin_ORBGRAND1(H,max_query,y_soft,even_code)

    % Hard demodulate
    y_demod = (y_soft<0);
    [k,n]=size(H);
    k=n-k;
    demod_parity = mod(sum(y_demod),2);

    n_guesses = 1;
    err_vec = zeros(1,n);

    %Abandon default values
    y_decoded = -1*ones(size(y_demod));

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

    % Calculate the syndrome
    Hy = mod(H*y_demod',2);

    % If the code is odd or code is even and it is an even error, the first
    % query is the demod strig
    if even_code==0 || (even_code && demod_parity==0)
        % Sum of the probability of guesses
        sPG=pG1;
        pG=pG1;
        % First query is demodulated string     
        if Hy==zeros(size(Hy))
            y_decoded = y_demod;
            abandoned = 0;
            if even_code == 0
                APP = pG/(pG+(1-sPG)*(2^k-1)/(2^n-n_guesses));
            else
                APP = pG/(pG+(1-sPG)*(2^k-1)/(2^(n-1)-n_guesses));
            end
            return;
        end
    end
    
    % If the demodulated string is not in the codebook, decode.
    [L, ind_order] = sort(abs_LLR,'ascend');
    % % Reorder the a priori bit flip probabilities
    % p_demod = p_demod(ind_order);

    % Inverse sort order
    inv_perm = zeros(1,length(ind_order));
    for ii=1:length(ind_order)
           inv_perm(ind_order(ii))=ii;    
    end

    % Slope
    beta = (L(round(n/2))-L(1))/(round(n/2)-1);
    % Intercept
    c = max(round(L(1)/beta-1),0);

    % This is the H columns reordered to put in ML order
    test_H = H(:,ind_order);
    % Total starting weight
    wt=c+1;
    while n_guesses<max_query && wt<=c*n+n*(n+1)/2
        % Hamming weight
        w=max(1,ceil((1+2*(n+c)-sqrt((1+2*(n+c))^2-8*wt))/2));
        % If the code is even and the parity is incorrect, increment 
        if even_code && rem(w,2)~=demod_parity
            w = w+1;
        end
        while w<=n
            % Logistic weight
            W = wt-c*w;
            if W<w*(w+1)/2
                break;
            else
                % Make error vectors
                % Internally converts W and n to W' and n'.
                noise_locations = landslide(W,w,n); 
                % If the code is odd and this is the first query

                
                % For each error vector
                for jj=1:size(noise_locations,1)
                    n_guesses = n_guesses +1;
                    err_vec = zeros(1,n);
                    indexes = noise_locations(jj,:);
                    err_vec(indexes)=1;
                    pG = pG1*prod(p_demod(ind_order(indexes))./(1-p_demod(ind_order(indexes))));
                    sPG = sPG+pG;

                    if (Hy == mod(test_H*err_vec',2))
                        err_vec = err_vec(inv_perm);
                        y_decoded = mod(y_demod-err_vec,2);
                        abandoned = 0;
                        if even_code == 0
                            APP = pG/(pG+(1-sPG)*(2^k-1)/(2^n-n_guesses));
                        else
                            APP = pG/(pG+(1-sPG)*(2^k-1)/(2^(n-1)-n_guesses));
                        end
                        return;
                    end
                end
            end
            % Increment Hamming weight 
            w=w+1;
            % If the code is even
            if even_code && rem(w,2)~=demod_parity
                w = w+1;
            end
        end
        wt=wt+1;
    end

    % If we max out on queries or total weight
    abandoned = 1;
    err_vec = zeros(size(y_demod));
end


% With W being the target logistic weight, w being the Hamming weight and n
% being the length of the string, W1 = W-w(w+1)/2 and n1 = n-w.

function z = landslide(W,w,n)

    W1=W-w*(w+1)/2;
    n1=n-w;
    % Create the first integer partition
    jj=1;
    % Start with empty vector and breaking at first index
    u = zeros(1,w);
    k=1;
    u = mountain_build(u,k,w,W1,n1);
    z(jj,:)=u;
    % Evaluate drops
    d=circshift(u,-1)-u;
    d(w)=0;
    % Evaluate accumuated drops
    D = cumsum(d,'reverse');
    % Each loop generates a new integer partition
    while D(1)>=2
        % Find the last index with an accumulated drop >=2
        k=find(D>=2,1,'last');
        % Increase its index by one.
        u(k)=u(k)+1;
        u = mountain_build(u,k,w,W1,n1);
        % Record the partition
        jj=jj+1;
        z(jj,:)=u;
        % Evaluate drops
        d=circshift(u,-1)-u;
        d(w)=0;
        % Evaluate acumuated drops
        D = cumsum(d,'reverse');
    end
    z = z + repmat([1:w],size(z,1),1);
end

function u = mountain_build(u,k,w,W1,n1)
    u(k+1:w) = u(k)*ones(1,w-k);
    W2 = W1-sum(u);
    q = floor(W2/(n1-u(k)));
    r = W2-q*(n1-u(k));
    if q ~= 0
        u(w-q+1:w)=n1*ones(1,q);
    end
    if w-q>0
    	u(w-q)=u(w-q)+r;
    end
end
