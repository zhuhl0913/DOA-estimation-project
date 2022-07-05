clear;
%% Load data





% load('Observation_wb.mat');                        % load data



% N1=length(X);

[X1,~]=audioread('data1/data1_01.wav');
[X2,~]=audioread('data1/data1_02.wav');
[X3,~]=audioread('data1/data1_03.wav');
[X4,fs]=audioread('data1/data1_04.wav');
[N1,b]=size(X1);

soundsc(X1,16000);


fra_per_win = 2123;                               %每个window中frame的数量
n_window = fix(N1/fra_per_win);                        %window个数
%每个window wlen个采样点
N=n_window*fra_per_win;                         %舍去多余的采样


X=zeros(N,4);
X(1:N,1)=X1(1:N,1);
X(1:N,2)=X2(1:N,1);
X(1:N,3)=X3(1:N,1);
X(1:N,4)=X4(1:N,1);

%% Plot waveform

wave_form = X;                                 % select one sensor's data as data sample

ds = 1/fs;                             
N = length(wave_form);
t = 0:ds:length(X)/fs-ds;                                    % set the x-axis
f = (-N/2:N/2-1)*fs/N;
Y = fftshift(fft(wave_form))/N;
subplot(2,1,1);plot(t,wave_form); xlabel('t(s)');ylabel('y(t)')
subplot(2,1,2);stem(f,abs(Y),'filled');axis([-8500 8500 0 inf]);
xlabel('f(Hz)');ylabel('Y(f)')
%[m,index_index] = max(abs(Y))
%% Array setup

[Frame,nSensors] = size(X);                            
J = nSensors;                                      % number of sensors
dx = 0.025;                                        % inter-sensor distance in x direction (m)
dy = 0;                                            % sensor distance in y direction (m)
c = 340;                                           % sound velocity  (m)
n_source = 2;                                      % number of sources
Index = linspace(0,J-1,J);
p = (-(J-1)/2 + Index.') * [dx dy];                % sensor position
%% Plot sensor positions

linspec = {'rx','MarkerSize',12,'LineWidth',2};
figure
plot(p(:,1),p(:,2),linspec{:});  
title('Sensor positions');
xlabel('x position in meters');
ylabel('y position in meters');
disp('The four microphones are ready !');
%% STFT（Short Time Fourier Transformation）



matrix_STFT = zeros([fra_per_win,n_window,nSensors]); %预分配一个 窗口长度*窗口数量*4 的零矩阵

for n = 1:nSensors                                      %将矩阵分割为 窗口长度*窗口数量*4 的矩阵
    for x = 1:n_window
    matrix_STFT(:,x,n) = X(fra_per_win*(x-1)+(1:fra_per_win),n);
    end
end


matrix_STFT_result = fft(matrix_STFT);
matrix_STFT_result = matrix_STFT_result(1:(fra_per_win/2+1),:,:);
%% Evaluate the general denominator

stride = 1;                                           % determine the angular resolution(deg)
theta = -90:stride:90;                                  % grid
adding_matrix = zeros(4,4);
denominator_temp = zeros(size(theta,2),1);
linear_frequency_distribution = linspace(0,fs,fra_per_win);
for n = 1:size(matrix_STFT_result,1)               
    temping_matrix = permute(matrix_STFT_result(n,:,:),[2 3 1]) ;
    transpose_matrix = temping_matrix.';
    adding_matrix = transpose_matrix*transpose_matrix';
        
    R_x = (1/n_window).*adding_matrix;                             % autocorrelation estimate
    v = [sin(theta*pi/180);-cos(theta*pi/180)];            % direction vector  
    a_theta = exp(-1j*2*pi*linear_frequency_distribution(n)*(p*v)./c);                  % steer vector
        
    % implement eigen-decomposition
    [V,D] = eig(R_x);
    for x = 1:J
        temp(1,x)=D(x,x);
    end
    V(J+1,:) = temp(1,:);
    V_temp = sortrows(V.',J+1);
    V_temp = V_temp.';
    V = V_temp(1:J,:);
    Un = V(:,1:J-n_source);
    denominator_temp = denominator_temp+diag(a_theta'*(Un*Un')*a_theta);
    
end
%% DoA estimation (MUSIC)                                    % noise subspace (columns are eigenvectors), size: J*(J-n_source)
P_sm = 1./denominator_temp;          % pseudo music power
%% Plot the MUSIC pseudo power spectrum
figure;
linspec = {'b-','LineWidth',2};
plot(theta, 10*log10(abs(P_sm)), linspec{:});
title('MUSIC pseudo power spectrum')
xlabel('Angle in [degrees]');
ylabel('Power spectrum in [dB]');
xlim([-90,90]);

%% Find the local maximum and visualization
P_middle = abs(P_sm(2:end-1));
P_front = abs(P_sm(1:end-2));
P_back = abs(P_sm(3:end));
logic_front = (P_middle - P_front)>0;
logic_back = (P_middle - P_back)>0;
logic = logic_front & logic_back;
P_middle(~logic) = min(P_middle);
P_local = [abs(P_sm(1));P_middle;abs(P_sm(end))];
[~,doa_Idx] = maxk(P_local,n_source);
doa = theta(doa_Idx);
[~,minIdx] = min(abs(doa));
doa_source = doa(minIdx);
[~,maxIdx] = max(abs(doa));
interfer = doa(maxIdx);
disp(['The desired source DOA with MUSIC is: ',num2str(doa_source),' deg']);
disp(['The interfering DOA with MUSIC is: ',num2str(interfer),' deg']);


