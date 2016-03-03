clear all
close all

%% specify a few parameters
%--------------------------------------------------------------------------

% number of code symbols
n=10000;

% channel parameter
%epsilon = 10.^[-1:-0.1:-2];
epsilon = 0.1:0.05:0.8;

% number of simulation runs
N_simulations = 100;


% two degree distributions
%lambda(x)
Distr_1A.degree   =  [ 3 ];
Distr_1A.fraction =  [ 1 ]; 
Distr_1A.fraction =  Distr_1A.fraction/sum(Distr_1A.fraction); %lambda_i
%rho(x)
Distr_2A.degree   =  [ 6 ];
Distr_2A.fraction =  [ 1 ];
Distr_2A.fraction =  Distr_2A.fraction/sum(Distr_2A.fraction); %rho_i

% number of iterations
N_iterations = 20;


% determine a few parameters in order to be able to set up the decoder...
%--------------------------------------------------------------------------

% Derive another two degree distributions
Distr_1B.degree   = Distr_1A.degree;
Distr_1B.fraction = Distr_1A.fraction./Distr_1A.degree; %fraction lambda_i/i
Distr_1B.fraction = Distr_1B.fraction/sum(Distr_1B.fraction); %lambda tilde

% Adjust edge distribution on variable nodes
N_node1 = round(Distr_1B.fraction*n); %array of number of variable nodes for each degree
if sum(N_node1)<n
    % add some variable nodes where there are the least
    [m,ind] = min(N_node1);
    N_node1(ind) = N_node1(ind) + (n-sum(N_node1));
elseif sum(N_node1)>n
    % remove some variable nodes where there are the most
    [m,ind] = max(N_node1);
    N_node1(ind) = N_node1(ind) - (sum(N_node1)-n);
end


% interleaver length, interleaver, and deinterleaver to distribute edges
% between variable-nodes and check-nodes
N             = N_node1 * Distr_1A.degree'; %Total number of edges going out of variable nodes
interleaver   = randperm(N);
deinterleaver = zeros(1,N);
for ii = 1:N
    deinterleaver(interleaver(ii)) = ii;
end


% Try to adjust edge distribution of check nodes
N_node2 = floor(Distr_2A.fraction./Distr_2A.degree*N);
degree_tmp =  N - N_node2*Distr_2A.degree'; %error of number of check nodes edges
if degree_tmp > 0 % if edges are missing
    if any(Distr_2A.degree==degree_tmp)
        % if adding 1 check node of a particular degree solve the issue, do it
        ind = find(Distr_2A.degree==degree_tmp,1,'first');
        N_node2(ind) = N_node2(ind) + 1;
    else
        % otherwise, create a new check-node degree equal to the number of missing edges
        Distr_2A.degree = [ Distr_2A.degree, degree_tmp];
        N_node2 = [ N_node2, 1];
    end
end

% Recalculate rho(x) after the adjustment
Distr_2B.degree   = Distr_2A.degree;
Distr_2B.fraction = N_node2/sum(N_node2); %rho_i



%% start the simulation
%--------------------------------------------------------------------------

% initialize BER
BER = zeros(1,length(epsilon));

tic
for ii_sim = 1:N_simulations
    for ii_e = 1:length(epsilon)
        
        % generate a codeword (BPSK)
        c = ones(1,n); %all-zero codeword
        
        % transmission over a noisy channel
        %e = 1-2*(rand(1,n)<epsilon(ii_e)); %BSC errors
        e = 1-(rand(1,n)<epsilon(ii_e)); %BEC erasures
        y = c.*e; %apply errors to codeword
                
        % ----------------------------------------
        % iterative decoding
        % ----------------------------------------
        
        % initialize the messages exchanged by the decoders
        in1 = zeros(1,N);
        in2 = zeros(1,N);
        out1 =zeros(1,N);
        out2 =zeros(1,N);
        
        
        for ii_it = 1:N_iterations
            % ----------------------------------------
            % run the decoders of type 1
            % ----------------------------------------
            
            pointer_A = 1; %input from channel pointer
            pointer_B = 1; %edges pointer
            
            for ii_d1 = 1:length(N_node1) %for each variable-node degree
                for ii_dec_d1 = 1:N_node1(ii_d1) %for each variable node
                    
                    % degree of the current decoder
                    degree1 = Distr_1B.degree(ii_d1);
                    
%                     % run the local decoder
%                     [ out1(pointer_B+[0:degree1-1]), c_hat(pointer_A) ]...
%                         = decoder_1( degree1, y(pointer_A), ...
%                         in1(pointer_B + [0:degree1-1]));
                    
                    % run the local decoder of BEC
                    [ out1(pointer_B+[0:degree1-1]), c_hat(pointer_A) ]...
                        = decoder_1belief( degree1, y(pointer_A), ...
                        in1(pointer_B + [0:degree1-1]));
                    
                    % increment "pointers"
                    pointer_A = pointer_A + 1;
                    pointer_B = pointer_B + degree1;
                    
                end
            end
            
            % interleave
            in2 = out1(interleaver);
            
            % ----------------------------------------
            % run the decoders of type 2
            % ----------------------------------------
            
            pointer_C = 1;
            for ii_d2 = 1:length(N_node2) %for each check-node degree 
                for ii_dec_d2 = 1:N_node2(ii_d2) %for each check-node
                    
                    % degree of the current decoder
                    degree2 = Distr_2B.degree(ii_d2);
                    
%                     % run the local decoder
%                     out2(pointer_C+[0:degree2-1]) = ...
%                         decoder_2( degree2, in2(pointer_C + [0:degree2-1]) );
                    
                    % run the local decoder
                    out2(pointer_C+[0:degree2-1]) = ...
                        decoder_2belief( degree2, in2(pointer_C + [0:degree2-1]) );
                    
                    % increment "pointer"
                    pointer_C = pointer_C + degree2;
                    
                end
            end
            
            % deinterleave
            in1 = out2(deinterleaver);
            
            % measure BER and stop if everything is correct
            ber = sum(c_hat~=1)/n;
            if ber == 0
                break
            end
        end
        
        % calculate instantaneous averaged BER
        BER(ii_e) = ( BER(ii_e)*(ii_sim-1) + ber )/ii_sim;
        
        semilogy(epsilon,BER)
        %axis([0,0.1,10^-4, 1])
        title([' Number of Simulations : ' num2str(ii_sim)])
        xlabel('epsilon')
        ylabel('BER')
        grid on
        drawnow
    end
end
T = toc;
disp([' ### Processing time : ' num2str(T) ' s'])


