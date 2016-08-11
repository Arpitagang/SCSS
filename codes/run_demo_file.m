clc, clear, close all

addpath(['istft']);

param.wlen = 512;                   % length of windowing signal
param.h = param.wlen/4;             % hop size
param.nfft = param.wlen;            % number of FFT points

% reading training data
[tr_clip{1}, param.fs]=audioread(['data',filesep,'train_file_female.wav']);
[tr_clip{2}, param.fs]=audioread(['data',filesep,'train_file_male.wav']);

% STFT of training data
Tr{1} = stft(tr_clip{1}, param.wlen, param.h, param.nfft, param.fs);
Tr{2} = stft(tr_clip{2}, param.wlen, param.h, param.nfft, param.fs);
Tr{1} = abs(Tr{1});
Tr{2} = abs(Tr{2});

% reading test data
[test{1}, ~]=audioread(['data',filesep,'test_file_female.wav']);
[test{2}, param.fs]=audioread(['data',filesep,'test_file_male.wav']);

% making signal to signal ratio 0dB
len=min(length(test{1}),length(test{2}));
test{1}=test{1}(1:len);
test{2}=test{2}(1:len);
test{2} = test{2}/std(test{2});
test{2} = test{2}*std(test{1});

x = test{1} + test{2};                  % mixed signal   

X = stft(x, param.wlen, param.h, param.nfft, param.fs);         % STFT of mixed signal

param.lambda = 100;  % regularisation factor

% structure of seeds for initialization
seeds.b1 = rng;
seeds.b2 = rng;
seeds.g1 = rng;
seeds.g2 = rng;
seeds.g = rng;

[SDR1,SIR1,SAR1,rec_test_no_fw] = sep_without_framework(Tr,test,X,param,seeds);
disp('*********performance without application of framework*********')
SDR = SDR1 
SIR = SIR1 
SAR = SAR1
disp('**************************************************************')

%% 

param.r_th = 3;               % threshold for error ratio

param.max_dim = 60;           % maximum dimension for binary search
param.min_dim = 20;           % minimum dimension for binary search

% separating one source at a time
% first source = source, second source = interferer
[SDR2(1,:),SIR2(1,:),SAR2(1,:),rec_test_fw{1}] = sep_with_framework(Tr{1},Tr{2},test{1},test{2},X,param,seeds);

% first source = source, second source = interferer
[SDR2(2,:),SIR2(2,:),SAR2(2,:),rec_test_fw{2}] = sep_with_framework(Tr{2},Tr{1},test{2},test{1},X,param,seeds);

disp('*********performance with application of framework*********')
SDR = SDR2
SIR = SIR2
SAR = SAR2
disp('***********************************************************')
