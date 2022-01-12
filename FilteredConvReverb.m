clear all
close all


[x, fs] = audioread ('drums.wav'); 
[h1] = audioread('ChurchImpulse48k.wav'); 
h1 = h1(:,1);


% z = zeros(length(x),1);
% z(1) = 1;
% 
% x = z;

% ------------------------------------------------------------------
% Prepare impulse signal
% ------------------------------------------------------------------

% Remove any silence from the start 
    framesize = 0.01*fs; % 0.01 per frame
    N = length(h1);
    nframes = floor(N/framesize);
    h_silenced = zeros(N,1);
    inx = 0;
    thresh = 0.001;
    
    for k=1:nframes
       frame = h1((k-1)*framesize+1 : framesize*k);
       maximum = max(frame);
        %only append signal at amplitude >0.02
       if(maximum > thresh)
            inx = inx+1;
            h_silenced((inx-1)*framesize+1 : framesize*inx) = frame;
       end
    end
    
% shorten impulse
h = h_silenced(1:fs);  

% ------------------------------------------------------------------
% Mixing Params
% -----------------------------------------------------------------
VA = 0.9; % verb amount
dry=0.5;
wet = 1-dry;

% ------------------------------------------------------------------
% Delay Line Params
% -----------------------------------------------------------------

D = length(h)/16;
DDL = zeros(1,D);
rwInx = 1;
fb = 0.9;
DelayDry = 0.9;
DelayWet = 1 - DelayDry;

% ------------------------------------------------------------------
% Cascaded Filters Params
% -----------------------------------------------------------------

LowGain = 0;
HighGain = 0;
LowCutoff = 1000;
HighCutoff = 15000;


% ------------------------------------------------------------------
% Filter Preperation & Design Equations
% -----------------------------------------------------------------

% Notmalise low params
GL = 10^(LowGain/20);            
wcL = 2*pi*LowCutoff/fs; 

% normalise high params
GH = 10^(HighGain/20);
wcH = 2*pi*HighCutoff/fs;

% Low shelve coefficients
b0Low = GL*tan(wcL/2) + sqrt(GL);
b1Low = GL*tan(wcL/2) - sqrt(GL);
a0Low = tan(wcL/2) + sqrt(GL);
a1Low = tan(wcL/2) - sqrt(GL);

% High shelve coefficients
b0High = sqrt(GH)*tan(wcH/2) + GH;
b1High = sqrt(GH)*tan(wcH/2) - GH;
a0High = sqrt(GH)*tan(wcH/2) + 1;
a1High = sqrt(GH)*tan(wcH/2) - 1;

% Normalising coeficcients
b0Low = b0Low/a0Low;
b1Low = b1Low/a0Low;
a1Low = a1Low/a0Low;
a0Low = a0Low/a0Low;

b0High = b0High/a0High;
b1High = b1High/a0High;
a1High = a1High/a0High;
a0High = a0High/a0High;

zX1L = 0;
zY1L = 0;

zX1H = 0;
zY1H = 0;

% ------------------------------------------------------------------
% Apply filter and delay Fx to impulse
% -----------------------------------------------------------------
for n=1:length(h)
    % Low shelve temp variable
    lowShelve = b0Low*h(n) + b1Low*zX1L - a1Low*zY1L;
    
    zX1L = h(n);
    zY1L = lowShelve;
    
   
    % High shelve temp variable
    highShelve = b0High*lowShelve + b1High*zX1H - a1High*zY1H; 
    
    zX1H = lowShelve;
    zY1H = highShelve;
    
    h_filt(n) = highShelve;
    % Read and write to delay line implementing feedback coefficient
    h_delay(n) = DelayDry*highShelve + DelayWet*DDL(rwInx);
    DDL(rwInx)= highShelve + fb * DDL(rwInx);
    
    % Implement circular buffer at length of delay line
if (rwInx>D)
    rwInx=1;
end
    
end

% Implement reverb amount. 
h_delay = h_delay*VA; 
h_delay = h_delay.';


% ------------------------------------------------------------------
% PRE-PROCESS 
% ------------------------------------------------------------------
l = 1024; % An arbitrary value to calulate the length of the FFT for new L sized chunks. 

ImpulseLength = length(h_delay); 
NFFT = 2^nextpow2(l+ImpulseLength-1); 
L = NFFT - ImpulseLength + 1; 

% ------------------------------------------------------------------
% Compute FFT of impulse response
% ------------------------------------------------------------------
h = [h_delay; zeros(NFFT - (ImpulseLength -1), 1)]; 
H = fft(h, NFFT); 

% ------------------------------------------------------------------
% Initiate Overlap & Output Buffer
% ------------------------------------------------------------------
OVERLAP = zeros((ImpulseLength - 1), 1); 
y = zeros(length(x) + ImpulseLength - 1, 1); 


% ------------------------------------------------------------------
% FFT Convolution Overlapp Add Loop
% ------------------------------------------------------------------
n = 0; 
while( n < (length(x) + ImpulseLength) / L) 
    
    % The end of the signal will not be consistent with other block sizes
    % so this is compensated for. 
     k = min(((n + 1)*(L)), length(x)); % returns the smaller size of the two numbers!!
     frameLength = (n * L) + 1 : k; 

     % Extract L sized chunk
     x_win = x(frameLength); 

     % Introduce overlap buffer before window 
     x_OLAP = [OVERLAP; x_win]; 
    
     % Pad the signal with zeros to even out signal lengths 
     x_zpad = [x_OLAP; zeros(length(H) - length(x_OLAP), 1)]; 

     OVERLAP = x_zpad(length(x_zpad) - (ImpulseLength - 2):length(x_zpad)); 

     % get DFT of the block 
     X_win = fft(x_zpad, NFFT); 
     
      % Window DFT
     Y_win = X_win .* H(:,1); 
     
      % Convolution 
     y_win = dry*real(ifft(X_win))+wet*real(ifft(Y_win)); 
     
     % Create y array with overlap add.
     y((n * L + 1: ((n + 1) * L) + 1), 1) = y_win(ImpulseLength-1 : length(y_win), 1); 
   
     n = n + 1; 
end 
%----------------------------------------------- 



% ------------------------------------------------------------------
% Normalise Output
% ------------------------------------------------------------------

y= y/max(max(abs(y))); 



