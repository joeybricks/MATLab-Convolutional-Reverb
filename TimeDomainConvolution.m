clear
close 
clc

tic 

[x,fs]= audioread('drums.wav'); %input wav file
[h,fs]= audioread('ChurchImpulse48k.wav');%impulse response

h = h(:,1); % make impulse mono
h = h (1:fs/2); % shorten impulse
N = (length(h)+length(x))-1;

% ------------------------------------------------------------------
% Prepare impulse signal
% ------------------------------------------------------------------

% Remove any silence from the start 
    frame_len = 0.01*fs; % 0.01 per frame
    N = length(h);
    num_frames = floor(N/frame_len);
    h_silenced = zeros(N,1);
    inx = 0;
    
    for k=1:num_frames
       frame = h((k-1)*frame_len+1 : frame_len*k);
       max_val = max(frame);
        %only append signal at amplitude >0.02
       if(max_val > 0.001)
            inx = inx+1;
            h_silenced((inx-1)*frame_len+1 : frame_len*inx) = frame;
       end
    end
    
    h = h_silenced;
    

% ------------------------------------------------------------------
% Zero Pad Impulse 
% ------------------------------------------------------------------
h = [h ;zeros(N-length(h),1)];
x = [x ;zeros(N-length(x),1)];

% % ------------------------------------------------------------------
% % Wet/Dry Mix Perams
% % ------------------------------------------------------------------
% dry = 0.5
% wet = 1-dry;

% % ------------------------------------------------------------------
% % Weighting h
% % ------------------------------------------------------------------
% VA = 0.7;

% ------------------------------------------------------------------
% Compute Convolution
% ------------------------------------------------------------------
convolvedSignal = conv(x,h);
% convolvedSignal = conv(x,VA*h);
% convolvedSignal = conv(dry*x,wet*h);
% convolvedSignal = conv(dry*x,wet*h*VA);
toc

plot(convolvedSignal,'k')
title('Convolutional with h = h(1:fs/2) Using conv() - Time domain')
hold on
plot(x,'r')
xlabel('Samples')
ylabel('Amplitude')
xlim([0 length(convolvedSignal)]);

