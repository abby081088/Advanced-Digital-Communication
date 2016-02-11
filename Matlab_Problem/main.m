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


  %------------------------------------------
  % Some parameters
  %------------------------------------------
  

  
  %------------------------------------------
  % parameters related to the data
  
  % number of data symbols
  N_symbols = 10;
  
  % number of data frames
  N_frames  = 1;
  
  % duration of one symbol
  T_sym    = 1;
  

  
  
  %------------------------------------------
  % parameters related to the channel model
    
  % resolution of the discrete-time implementation
  %  (note that all sampling times etc. should 
  %   be multiples of delta_t)
  delta_t  = 0.05;
  
  
  % SNR range (in dB)
  EbN0_dB = -10:2.5:25;
  
  
  % specify the ISI channel 
  Channel.Gains  = [ 2 -3/4 1i];
  Channel.Delays = [ 0.5 2 2.25];
  
  
  
  
  %------------------------------------------
  % parameters related to the receiver
  
  % sampling factor and sampling time
  m        = 2;
  T_sample = T_sym/m;
  
  
  
  
  %------------------------------------------
  % Define some filters
  %------------------------------------------
  
  
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
  
% $$$ T_F4     = 1;
% $$$ Filter_4 = [0; ones(round(T_F4/delta_t),1,1)];


  % generate filter 5
  T_F5     = T_F3 + T_F4;
  Filter_5 = conv(Filter_3,Filter_4)*delta_t;
 

  
  
  
  %------------------------------------------
  % Start the simulation
  %------------------------------------------
  
  
  for ii_F = 1:N_frames
    for ii_SNR = 1:length(EbN0_dB)
      

      
      %------------------------------------------
      % Step 1
      %------------------------------------------
      

      % generate vector 1
      Vec_1    = (2*(rand(N_symbols,2)>0.5) -1)*[1;1i];
      P        = 2;
      
      
      % generate signal 1
      Signal_1 = zeros((N_symbols-1)*T_sym/delta_t+1,1);
      Signal_1(1:round(T_sym/delta_t):end) = Vec_1/delta_t;
    

      % generate signal 2
      Signal_2 = conv(Signal_1,Filter_1)*delta_t;
  


      %------------------------------------------
      % Step 2
      %------------------------------------------
      

      % generate  signal 3
      Signal_3 = conv(Signal_2,Filter_2)*delta_t;
      

      % generate another signal
      % (note that the factor 1/sqrt(delta_t) is 
      %  needed to get the noise variance after 
      %  the receiver filter correct.)
      Eb = P/2 * (Filter_3'*Filter_3)*delta_t;
      g  = sqrt(Eb/10^(EbN0_dB(ii_SNR)/10));
      Signal_4 = g*1/sqrt(2)*(randn(length(Signal_3),2)*[1;1i])/sqrt(delta_t);
      
      
      % generate signal 5
      Signal_5 = Signal_3 + Signal_4;
      



      %------------------------------------------
      % Step 3
      %------------------------------------------
      

      % generate signal 6
      Signal_6 = conv(Filter_4,Signal_5)*delta_t;


      % generate vector 2  
      Offset = rem(T_F3,T_sym);
      Vec_2  = Signal_6([1+round(Offset/delta_t):round(T_sample/delta_t):end]);
             
  
      % generate a discrete filter
      d_filter = Filter_5([1+round(Offset/delta_t):round(T_sample/delta_t):end]);
      
      
      %------------------------------------------
      % some preprocessing
      %------------------------------------------
      
      % remove leading zeros (decision delay)
      i1 = find(d_filter ~= 0,1,'first');
      i2 = find(d_filter ~= 0,1,'last');
      
      d_filter = d_filter(i1:i2);
      Vec_2      = Vec_2(i1:end);
      




      %------------------------------------------
      % Step 4
      %------------------------------------------
      

      % observation interval
      L_o = length(d_filter);
      

      % zero-padding if observation interval is langer than the pulse
      if L_o>length(Vec_2)
	Vec_2 = [ Vec_2; zeros(L_o-length(Vec_2),1) ];
      end
      
      
      
      %------------------------------------------
      % Equalization
      %------------------------------------------
      

      %
      %  Here, you should implement the equalizers
      %  and measure the BER. The function GenerateMatrix.m 
      %  may be helpful for you.     
      %     
      
      
      
      
      
    end
  end
  
  
  %------------------------------------------
  % Plot results
  %------------------------------------------
  
  
  %
  %  Here, you should plot your results.
  %  Use semilogy(EbN0_dB,BER) for the plots 
  %  and do not forget the legend!   
  
  
  
  
  %------------------------------------------
  % Look at the waveforms
  %------------------------------------------
  figure(2)
  subplot(511), plot([0:length(Signal_2)-1]*delta_t,Signal_2),grid on
  x_max =  (length(Signal_2)-1)*delta_t;
  y_max = 1.2*max(abs(Signal_2));
  axis([0 x_max -y_max y_max ])
  
  subplot(512), plot([0:length(Signal_3)-1]*delta_t,Signal_3),grid on
  x_max =  (length(Signal_3)-1)*delta_t;
  y_max = 1.2*max(abs(Signal_3));
  axis([0 x_max -y_max y_max ])
  
  subplot(513), plot([0:length(Signal_5)-1]*delta_t,Signal_5),grid on
  x_max =  (length(Signal_5)-1)*delta_t;
  y_max = 1.2*max(abs(Signal_5));
  axis([0 x_max -y_max y_max ])
  
  subplot(514), plot([0:length(Signal_6)-1]*delta_t,Signal_6),grid on
  x_max =  (length(Signal_6)-1)*delta_t;
  y_max = 1.2*max(abs(Signal_6));
  axis([0 x_max -y_max y_max ])
  
  hold on
  stem(Offset + (i1-1+[0:length(Vec_2)-1])*T_sample,Vec_2 )
  hold off
  
  subplot(515), plot([0:length(Filter_5)-1]*delta_t,Filter_5),grid on
  
  
  %------------------------------------------
  %------------------------------------------
  
