clear all
close all
clc;


%% Radar Specifications
% Part of Uacitity Sensor Fusion Project
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Frequency of operation = 77GHz
% Max Range = 200m
% Range Resolution = 1 m
% Max Velocity = 100 m/s
%%%%%%%%%%%%%%%%%%%%%%%%%%%
Res = 1; % Resolution in m
Rmax = 200; % max range in m

%speed of light = 3e8
c = 3e8; % speed of light in m/s

%% Initial values - User Defined Range and Velocity of target
R = 110; % target position in m
v = -20; % target speed in m/s


%% FMCW Waveform Generation

% Design the FMCW waveform by giving the specs of each of its parameters.
% Calculate the Bandwidth (B), Chirp Time (Tchirp) and Slope (slope) of the FMCW
% chirp using the requirements above.

Bandwith = c / (2 * Res);       % Bandwith in 1/s or Hz
Tchirp = 5.5 * (2 * Rmax) / c;  % chirp time assuming B/Fs = 5.5 where Fs is ADC sampling Freq
Slope = Bandwith / Tchirp;      % slope of chirp frequency change in Hz/sec


%Operating carrier frequency of Radar 
fc= 77e9;             %carrier freq

                                                          
%The number of chirps in one sequence. Its ideal to have 2^ value for the ease of running the FFT
%for Doppler Estimation. 
Nd=128;                   % #of doppler cells OR #of sent periods % number of chirps

%The number of samples on each chirp. 
Nr=1024;                  %for length of time OR # of range cells

% Timestamp for running the displacement scenario for every sample on each
% chirp
t=linspace(0,Nd*Tchirp,Nr*Nd); %total time for samples


%Creating the vectors for Tx, Rx and Mix based on the total samples input.
Tx=zeros(1,length(t));      % transmitted signal
Rx=zeros(1,length(t));      % received signal
Mix = zeros(1,length(t));   % beat signal

%Similar vectors for range_covered and time delay.
r_t=zeros(1,length(t));
td=zeros(1,length(t));

%% Signal generation and Moving Target simulation
% Simulating the radar scenario over the time. 

for i=1:length(t)         
    %For each time stamp update the Range of the Target for constant velocity. 
    r_t(i) = R + v * t(i);
    td(i)  = 2 * r_t(i) / c;    % sec : ToF
    
    %For each time sample we need update the transmitted and received signal. 
    Tx(i)   = cos(2 * pi * (fc * t(i) + 0.5 * Slope * t(i) ^ 2));                   % a.u. : transmitted signal
    Rx(i)   = cos(2 * pi * (fc * (t(i) - td(i)) + 0.5 * Slope * (t(i) - td(i))^2)); % a.u. : received signal - Amplitude ignored
    
    %Now by mixing the Transmit and Receive generate the beat signal
    %This is done by element wise matrix multiplication of Transmit and
    %Receiver Signal
    Mix(i) = Tx(i) * Rx(i);
end


%% RANGE MEASUREMENT

%reshape the vector into Nr*Nd array. Nr and Nd here would also define the size of
%Range and Doppler FFT respectively.
Mix = reshape(Mix, [Nr, Nd]);

%run the FFT on the beat signal along the range bins dimension (Nr) and
%normalize.
fft_mix = fft(Mix, Nr, 1);

% Take the absolute value of FFT output
fft_mix = abs(fft_mix);
fft_mix = fft_mix(1:Nr/2+1);


%plotting the range
figure ('Name','Range from First FFT')
subplot(2,1,1)

% plot FFT output 
plot(fft_mix);
axis ([0 200 0 1]);


%% RANGE DOPPLER RESPONSE

% Range Doppler Map Generation.
% The output of the 2D FFT is an image that has reponse in the range and
% doppler FFT bins. So, it is important to convert the axis from bin sizes
% to range and doppler based on their Max values.

Mix=reshape(Mix,[Nr,Nd]);   % again here for completeness

% 2D FFT using the FFT size for both dimensions.
sig_fft2 = fft2(Mix,Nr,Nd);

% Taking just one side of signal from Range dimension.
sig_fft2 = sig_fft2(1:Nr/2,1:Nd);
sig_fft2 = fftshift (sig_fft2);
RDM = abs(sig_fft2);
RDM = 10*log10(RDM) ;

%use the surf function to plot the output of 2DFFT and to show axis in both dimensions
doppler_axis = linspace(-100,100,Nd);
range_axis = linspace(-200,200,Nr/2)*((Nr/2)/400);
figure,surf(doppler_axis,range_axis,RDM);


%% CFAR implementation

%Slide Window through the complete Range Doppler Map

%Select the number of Training Cells in both the dimensions.
Tr = 10; 
Td = 8;

%Select the number of Guard Cells in both dimensions around the Cell under test (CUT) for accurate estimation
Gr = 4;
Gd = 4;

% offset the threshold by SNR value in dB
offset = 6;

%Create a vector to store noise_level for each iteration on training cells
noise_level = zeros(1,1);


%iterative over the entire 2D fft map and compare the SNR of the training(neighborhood/local) cells around each cell to it

CFAR = RDM;     %create a copy
Number_of_cells = 4 * (Td + Gd + 1) * (Td + Gd + 1) - (Gr * Gd) - 1;    %from the square with edge size of (Td + Gd + 1) subtract the number of guard cells and the cell itself
for i = Tr + Gr + 1 : Nr/2 - (Gr + Tr)          % start from Tr+Gr+1 from edge and stop Tr + Gr + 1 from the other side edge
    for j = Td + Gd + 1 : Nd - (Gd + Td)        % same as above for the other direction
        
        local_noise = zeros(1,1);
        % Use RDM[x,y] as the matrix from the output of 2D FFT for implementing CFAR
        for p = i - (Tr + Gr) : i + Tr + Gr % Guaranteed that we dont go out of the 2D matrix because of above
            for q = j - (Td + Gd) : j + Td + Gd
                if (abs(i - p) > Gr || abs(j - q) > Gd)
                    local_noise = local_noise + db2pow(RDM(p, q));
                end
            end
        end
        threshold = pow2db(local_noise / Number_of_cells) + offset;
        
        if(RDM(i, j) < threshold)
            CFAR(i, j) = 0;
        else
            CFAR(i, j) = 1;
        end

    end
end

% To keep the map size same set the unvisited cells to 0. 
for i = 1 : Nr / 2
    for j = 1 : Nd
        if (i > (Tr + Gr)) & (j > (Td + Gd)) & (i < (Nr/ 2 - (Tr + Gd))) & (j < (Nd - (Td + Gd)))
            continue;
        end
        CFAR(i, j) = 0;
    end
end

%display the CFAR output for Range and Doppler Response output.
figure,surf(doppler_axis,range_axis,CFAR);
colorbar;
 
 