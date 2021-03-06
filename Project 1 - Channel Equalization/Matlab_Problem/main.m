%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                             %
%  Homework Project No.1                      %    
%  Advanced Digital Communications (EQ2410)   %
%  KTH/EES, Stockholm, Sweden                 %  
%  Period 3, 2013/14                          %
%                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  clear
  clf

  
  %% Some parameters
  %------------------------------------------------------------------------
  % parameters related to the data:
  
  % number of data symbols
  N_symbols = 500;
  
  % number of data frames
  N_frames  = 10;
  
  % duration of one symbol
  T_sym    = 1;
   
  %------------------------------------------------------------------------
  % parameters related to the channel model:
    
  % resolution of the discrete-time implementation
  %  (note that all sampling times etc. should 
  %   be multiples of delta_t)
  delta_t  = 0.05;
  
  % SNR range (in dB)
  EbN0_dB = -10:2.5:25;
  %EbN0_dB = 50; % only one SNR for debug
  
  % specify the ISI channel 
  Channel.Gains  = [ 2 -3/4 1i];
  Channel.Delays = [ 0.5 2 2.25];

  %------------------------------------------------------------------------
  % parameters related to the receiver:
  
  % sampling factor and sampling time
  m        = 2; 
  T_sample = T_sym/m; %m=1,2,4 to have T_sample multiple of delta_t 
  
  
  
  %% Define some filters
  %------------------------------------------------------------------------
  
  % generate filter 1
  T_F1     = T_sym;
  Filter_1 = ones(round(T_F1/delta_t),1,1);
  
  
  % generate filter 2
  T_F2     = max(Channel.Delays);
  Filter_2 = zeros(round(T_F2/delta_t)+1,1,1);
  Filter_2(Channel.Delays/delta_t+1) = Channel.Gains/delta_t;
  

  % generate filter 3
  T_F3     = T_F1 + T_F2;
  Filter_3 = conv(Filter_1,Filter_2)*delta_t;
  

  % generate filter 4 
  T_F4     = T_F3;
  Filter_4 = [0; conj(flipud(Filter_3))];
  
  %T_F4     = 1;
  %Filter_4 = [0; ones(round(T_F4/delta_t),1,1)];


  % generate filter 5
  T_F5     = T_F3 + T_F4;
  Filter_5 = conv(Filter_3,Filter_4)*delta_t;
 

  %% Start the simulation
  %------------------------------------------------------------------------
  
  BER_zf = zeros(length(EbN0_dB),1);
  BER_mmse = zeros(length(EbN0_dB),1);
  BER_dfe_zf = zeros(length(EbN0_dB),1);
  BER_dfe_mmse = zeros(length(EbN0_dB),1);
 
  for ii_F = 1:N_frames
    % Signals independant of SNR
    %----------------------------------------------------------------------
    % Step 1
    %----------------------------------------------------------------------
    % generate vector 1
    Vec_1    = (2*(rand(N_symbols,2)>0.5) -1)*[1;1i];
    P        = 2;      
      
    % generate signal 1
    Signal_1 = zeros((N_symbols-1)*T_sym/delta_t+1,1);
    Signal_1(1:round(T_sym/delta_t):end) = Vec_1/delta_t;

    % generate signal 2
    Signal_2 = conv(Signal_1,Filter_1)*delta_t;
    
    %----------------------------------------------------------------------
    % Step 2
    %----------------------------------------------------------------------
    % generate  signal 3
    Signal_3 = conv(Signal_2,Filter_2)*delta_t;

    for ii_SNR = 1:length(EbN0_dB)
      % generate another signal
      % (note that the factor 1/sqrt(delta_t) is 
      %  needed to get the noise variance after 
      %  the receiver filter correct.)
      Eb = P/2 * (Filter_3'*Filter_3)*delta_t;
      g  = sqrt(Eb/10^(EbN0_dB(ii_SNR)/10));
      Signal_4 = g*1/sqrt(2)*(randn(length(Signal_3),2)*[1;1i])/sqrt(delta_t);
      
      
      % generate signal 5
      Signal_5 = Signal_3 + Signal_4;
      
      %--------------------------------------------------------------------
      % Step 3
      %--------------------------------------------------------------------
      % generate signal 6
      Signal_6 = conv(Filter_4,Signal_5)*delta_t;

      % generate vector 2  
      Offset = rem(T_F3,T_sym); 
      %Remove small time at the begining of Signal_6 so that the duration of
      %Signal_6 is an integer multiple of symbol duration
      Vec_2  = Signal_6([1+round(Offset/delta_t):round(T_sample/delta_t):end]);
  
      % generate a discrete filter
      d_filter = Filter_5([1+round(Offset/delta_t):round(T_sample/delta_t):end]);
      
      %--------------------------------------------------------------------
      % some preprocessing
      %--------------------------------------------------------------------
      
      % remove leading zeros (decision delay)
      i1 = find(d_filter ~= 0,1,'first');
      i2 = find(d_filter ~= 0,1,'last');
      
      d_filter = d_filter(i1:i2);
      Vec_2      = Vec_2(i1:end);

      %--------------------------------------------------------------------
      % Step 4
      %--------------------------------------------------------------------
      

      % observation interval (i2-i1+1)
      L_o = length(d_filter);
      

      % zero-padding if observation interval is langer than the pulse
      if L_o>length(Vec_2)
        Vec_2 = [ Vec_2; zeros(L_o-length(Vec_2),1) ];
      end
      
      
      
      %--------------------------------------------------------------------
      % Equalization
      %--------------------------------------------------------------------
      Constellation = [1+1i,1-1i,-1-1i,-1+1i];
      k0=1;
      k1=floor((L_o-1)/m);
      k2=k1;
      U = GenerateMatrix(d_filter,L_o,k0,k1,k2,m);
      e = (1:k1+k2+1==(k1+1))';
      
      %--------------------------------------------------------------------
      % Zero-forcing equalizer
      if m >1 %ZF doesn't exist for m=1
          c_zf = (U/(U'*U))*e;
      
          
          %One can test this correlator with:
          %k=2;c_zf'*[zeros(k*m,1);d_filter(1:length(d_filter)-k*m)]
          %or k=3;c_zf'*[d_Zfilter(k*m+1:length(d_filter));zeros(k*m,1)]

          errors = 0;
          Symbols_zf = zeros(N_symbols,1);
          for ii_n = 0:N_symbols-1
              %Decision variable
              Z = c_zf'*Vec_2(1+ii_n*m:L_o+ii_n*m);

              %Hard decision
              dist = abs(Constellation - Z);
              [~,hard_dec] = min(dist);
              Symbols_zf(1+ii_n) = Constellation(hard_dec);

              % Calculate BER (we assume Gray coding is used)
              if(abs(Symbols_zf(1+ii_n)-Vec_1(1+ii_n))==2)
                  errors=errors+1;
              elseif(abs(Symbols_zf(1+ii_n)-Vec_1(1+ii_n))>2)
                  errors=errors+2;
              end
          end
          BER_zf(ii_SNR) =  BER_zf(ii_SNR) + errors/(2*length(Vec_1));
      end
      
      %--------------------------------------------------------------------
      % Linear MMSE equalizer
      P_s = 2;
      Sigma2noise = Eb/10^(EbN0_dB(ii_SNR)/10);
      autocorr_filter4 = mean(abs(Filter_4).^2)*autocorr([Filter_4;zeros(L_o*round(T_sample/delta_t),1)],L_o*round(T_sample/delta_t)-1);
      C_w=toeplitz(2*Sigma2noise*autocorr_filter4(1:round(T_sample/delta_t):end));
      
      R = P_s*(U*U')+C_w;
      p = P_s^2*U*e;
      c_mmse=R\p;
      
      errors = 0;
      Symbols_mmse = zeros(N_symbols,1);
      for ii_n = 0:N_symbols-1
          %Decision variable
          Z = c_mmse'*Vec_2(1+ii_n*m:L_o+ii_n*m);
          
          %Hard decision
          dist = abs(Constellation - Z);
          [~,hard_dec] = min(dist);
          Symbols_mmse(1+ii_n) = Constellation(hard_dec);
          
          % Calculate BER (we assume Gray coding is used)
          if(abs(Symbols_mmse(1+ii_n)-Vec_1(1+ii_n))==2)
              errors=errors+1;
          elseif(abs(Symbols_mmse(1+ii_n)-Vec_1(1+ii_n))>2)
              errors=errors+2;
          end
      end
      BER_mmse(ii_SNR) =  BER_mmse(ii_SNR) + errors/(2*length(Vec_1));

      %--------------------------------------------------------------------
      % Decision feedback equalization

      U_dfe = U(:,k1+1:end);
      e_dfe = (1:k2+1==1)';
      
      % Zero-forcing equalizer
      c_dfe_zf_ff = (U_dfe/(U_dfe'*U_dfe))*e_dfe; %feedforward correlator
      c_dfe_zf_fb = -c_dfe_zf_ff'*U(:,1:k1);
      
      errors = 0;
      Symbols_dfe_zf = zeros(N_symbols,1);
      for ii_n = 0:N_symbols-1
          %Decision variable
          Symbols_dfe_zf_padded=[zeros(k1,1);Symbols_dfe_zf]; %add leading zeros
          Z = c_dfe_zf_fb*Symbols_dfe_zf_padded(ii_n+1:ii_n+length(c_dfe_zf_fb),:) + c_dfe_zf_ff'*Vec_2(1+ii_n*m:L_o+ii_n*m);

          %Hard decision
          dist = abs(Constellation - Z);
          [~,hard_dec] = min(dist);
          Symbols_dfe_zf(1+ii_n) = Constellation(hard_dec);

          % Calculate BER (we assume Gray coding is used)
          if(abs(Symbols_dfe_zf(1+ii_n)-Vec_1(1+ii_n))==2)
              errors=errors+1;
          elseif(abs(Symbols_dfe_zf(1+ii_n)-Vec_1(1+ii_n))>2)
              errors=errors+2;
          end
      end
      BER_dfe_zf(ii_SNR) =  BER_dfe_zf(ii_SNR) + errors/(2*length(Vec_1));
    
      
      % MMSE
      R_dfe = P_s*(U_dfe*U_dfe')+C_w;
      p_dfe = P_s^2*U_dfe*e_dfe;
      c_dfe_mmse_ff=R_dfe\p_dfe;
      c_dfe_mmse_fb = -c_dfe_mmse_ff'*U(:,1:k1);
      
      errors = 0;
      Symbols_dfe_mmse = zeros(N_symbols,1);
      for ii_n = 0:N_symbols-1
          %Decision variable
          Symbols_dfe_mmse_padded=[zeros(k1,1);Symbols_dfe_mmse]; %add leading zeros
          Z = c_dfe_mmse_fb*Symbols_dfe_mmse_padded(ii_n+1:ii_n+length(c_dfe_mmse_fb),:) + c_dfe_mmse_ff'*Vec_2(1+ii_n*m:L_o+ii_n*m);

          %Hard decision
          dist = abs(Constellation - Z);
          [~,hard_dec] = min(dist);
          Symbols_dfe_mmse(1+ii_n) = Constellation(hard_dec);
          
          % Calculate BER (we assume Gray coding is used)
          if(abs(Symbols_dfe_mmse(1+ii_n)-Vec_1(1+ii_n))==2)
              errors=errors+1;
          elseif(abs(Symbols_dfe_mmse(1+ii_n)-Vec_1(1+ii_n))>2)
              errors=errors+2;
          end
      end
      BER_dfe_mmse(ii_SNR) =  BER_dfe_mmse(ii_SNR) + errors/(2*length(Vec_1));
      
    end
  end
  
  BER_zf = BER_zf/N_frames;
  BER_mmse = BER_mmse/N_frames;
  BER_dfe_zf = BER_dfe_zf/N_frames;
  BER_dfe_mmse = BER_dfe_mmse/N_frames;
  
  %------------------------------------------------------------------------
  %% Plot results
  %
  %  Here, you should plot your results.
  %  Use semilogy(EbN0_dB,BER) for the plots 
  %  and do not forget the legend!   
  
  figure(1);
  semilogy(EbN0_dB,BER_zf)
  hold on
  semilogy(EbN0_dB,BER_mmse)
  semilogy(EbN0_dB,BER_dfe_zf)
  semilogy(EbN0_dB,BER_dfe_mmse)
  hold off;
  xlabel('SNR (dB)')
  ylabel('BER (log scale)')
  title('BER of egalization algorithms');
  legend('ZF','MMSE','ZF-DFE','MMSE-DFE');
  
  %------------------------------------------------------------------------
  %% Look at the waveforms
  figure(2)
  % Transmitted signal (real part)
  subplot(511), plot([0:length(Signal_2)-1]*delta_t,Signal_2),grid on
  x_max =  (length(Signal_2)-1)*delta_t;
  y_max = 1.2*max(abs(Signal_2));
  axis([0 x_max -y_max y_max ])
  
  % Signal passed through the channel (real part)
  subplot(512), plot([0:length(Signal_3)-1]*delta_t,Signal_3),grid on
  x_max =  (length(Signal_3)-1)*delta_t;
  y_max = 1.2*max(abs(Signal_3));
  axis([0 x_max -y_max y_max ])
  
  % Signal passed through the channel + noise (real part)
  subplot(513), plot([0:length(Signal_5)-1]*delta_t,Signal_5),grid on
  x_max =  (length(Signal_5)-1)*delta_t;
  y_max = 1.2*max(abs(Signal_5));
  axis([0 x_max -y_max y_max ])
  
  % Signal after the receiver filter (real part)
  subplot(514), plot([0:length(Signal_6)-1]*delta_t,Signal_6),grid on
  x_max =  (length(Signal_6)-1)*delta_t;
  y_max = 1.2*max(abs(Signal_6));
  axis([0 x_max -y_max y_max ])
  hold on
  % Plot circles representing sample time of signal (real part)
  stem(Offset + (i1-1+[0:length(Vec_2)-1])*T_sample,Vec_2 )
  hold off
   
  % Overall system response (transmitter*channel*receiver filter) (real part)
  subplot(515), plot([0:length(Filter_5)-1]*delta_t,Filter_5),grid on
  
  
  %------------------------------------------------------------------------
  %% Plot the filters gain
  figure(3)  
  % Transmitter filter
  subplot(511), plot([0:length(Filter_1)-1]*delta_t,abs(Filter_1)),grid on
  x_max =  (length(Filter_1)-1)*delta_t;
  y_max = 1.2*max(abs(Filter_1));
  axis([0 x_max -y_max y_max ])
  
  % Channel filter
  subplot(512), plot([0:length(Filter_2)-1]*delta_t,abs(Filter_2)),grid on
  x_max =  (length(Filter_2)-1)*delta_t;
  y_max = 1.2*max(abs(Filter_2));
  axis([0 x_max -y_max y_max ])
  
  % Tansmitter*channel
  subplot(513), plot([0:length(Filter_3)-1]*delta_t,abs(Filter_3)),grid on
  x_max =  (length(Filter_3)-1)*delta_t;
  y_max = 1.2*max(abs(Filter_3));
  axis([0 x_max -y_max y_max ])
  
  % Receiver filter (matched filter)
  subplot(514), plot([0:length(Filter_4)-1]*delta_t,abs(Filter_4)),grid on
  x_max =  (length(Filter_4)-1)*delta_t;
  y_max = 1.2*max(abs(Filter_4));
  axis([0 x_max -y_max y_max ])
  
  % Overall system response
  subplot(515), plot([0:length(Filter_5)-1]*delta_t,abs(Filter_5)),grid on
  x_max =  (length(Filter_5)-1)*delta_t;
  y_max = 1.2*max(abs(Filter_5));
  axis([0 x_max -y_max y_max ])
  hold on
  % Plot circles representing sample time of filter 
  stem(Offset + (i1-1+[0:length(d_filter)-1])*T_sample,abs(d_filter) )
  hold off
 
  %------------------------------------------------------------------------
  
