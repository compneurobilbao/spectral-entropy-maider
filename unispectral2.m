%% clear all
clc

tic
%% Unispectral analysis
Fs = 500;  % Sampling frequency  
temp_win_sz=500; % Temporal window size (temp_win_sz/Fs = seconds of )
temp_win_step_sz=1; %Number of moving samples between one temporal window and the next one

freq = 1:Fs/temp_win_sz:Fs/2;

freq_win_sz=5; % Entropies calc window size
  

% Load data
filePattern=fullfile('/path/to/data/*.txt');
files=dir(filePattern);
num_users=length(files);

% Initialize variables
ppgpain=zeros(num_users,length(freq)-freq_win_sz);
ppgnopain=zeros(num_users,length(freq)-freq_win_sz);
eeg1pain=zeros(num_users,length(freq)-freq_win_sz);
eeg1nopain=zeros(num_users,length(freq)-freq_win_sz);
eeg2pain=zeros(num_users,length(freq)-freq_win_sz);
eeg2nopain=zeros(num_users,length(freq)-freq_win_sz);


% for each user
 for u=1:num_users 
     filename=files(u);
     name=filename.name;
     X=load(name);
     PPG=X(:,1);
     EEG1=X(:,3);
     EEG2=X(:,4);
     % GSR=X(:,2);
     PULS=X(:,5);
       
    % Pain / No Pain area limits
    for n=1:length(PPG)
     if PULS(n)==0
            break;
     end
    end

    n_startpain=n+1000;
    for m=n_startpain:length(PPG)
        if PULS(m)==0
            break;
        end
    end

    m_startpain=m+1000;
    
   % Signals to analyze
    PPGNoPain=PPG(m_startpain:end,:);
    PPGPain=PPG(n_startpain:m,:);
    EEG1Pain=EEG1(n_startpain:m,:);    
    EEG1NoPain=EEG1(m_startpain:end,:);
    EEG2Pain=EEG2(n_startpain:m,:);     
    EEG2NoPain=EEG2(m_startpain:end,:);
        
    % Entropies calc
    for k=1:temp_win_step_sz:length(PPGPain)-temp_win_sz
        x=PPGPain(k:k+temp_win_sz);
        x=x.*blackman(length(x));
        
        % Compute power spectrum and normalize to probability density
        % function
        N = length(x);
        xdft = fft(x);
        xdft = xdft(1:(N/2+1));
        psdx = abs(xdft).^2;
        pdfx= psdx./sum(psdx) ; % Normalize to get a probability density function
        
        
        % Entropy Calculation
        for j = 1:length(freq)-freq_win_sz
           aux1=pdfx(j:j+freq_win_sz);
           entropy_ppgpain(k,j)=-sum(aux1.*log2(aux1 + 1e-12))/log2(length(aux1));
        end
    end
    
      
    for k=1:temp_win_step_sz:length(PPGNoPain)-temp_win_sz
        x=PPGNoPain(k:k+temp_win_sz);
        x=x.*blackman(length(x));
        % Compute power spectrum and normalize to probability density
        % function
        N = length(x);
        xdft = fft(x);
        xdft = xdft(1:(N/2+1));
        psdx = abs(xdft).^2;
        pdfx= psdx./sum(psdx) ; % Normalize to get a probability density function
        
        % Entropy Calculation
        for j = 1:length(freq)-freq_win_sz
           aux1=pdfx(j:j+freq_win_sz);
           entropy_ppgnopain(k,j)=-sum(aux1.*log2(aux1 + 1e-12))/log2(length(aux1));
        end
        
    end
    
    
     for k=1:temp_win_step_sz:length(EEG1Pain)-temp_win_sz
        x=EEG1Pain(k:k+temp_win_sz);
        x=x.*blackman(length(x));
        
        % Compute power spectrum and normalize to probability density
        % function
        N = length(x);
        xdft = fft(x);
        xdft = xdft(1:(N/2+1));
        psdx = abs(xdft).^2;
        pdfx= psdx./sum(psdx) ; % Normalize to get a probability density function
        
        % Entropy Calculation
        for j = 1:length(freq)-freq_win_sz
           aux1=pdfx(j:j+freq_win_sz);
           entropy_eeg1pain(k,j)=-sum(aux1.*log2(aux1 + 1e-12))/log2(length(aux1));
        end
        
     end
    
     for k=1:temp_win_step_sz:length(EEG1NoPain)-temp_win_sz
        x=EEG1NoPain(k:k+temp_win_sz);
        x=x.*blackman(length(x));
        
        % Compute power spectrum and normalize to probability density
        % function
        N = length(x);
        xdft = fft(x);
        xdft = xdft(1:(N/2+1));
        psdx = abs(xdft).^2;
        pdfx= psdx./sum(psdx) ; % Normalize to get a probability density function
        
        % Entropy Calculation
        for j = 1:length(freq)-freq_win_sz
           aux1=pdfx(j:j+freq_win_sz);
           entropy_eeg1nopain(k,j)=-sum(aux1.*log2(aux1 + 1e-12))/log2(length(aux1));
        end
        
     end
    
     for k=1:temp_win_step_sz:length(EEG2Pain)-temp_win_sz
        x=EEG2Pain(k:k+temp_win_sz);
        x=x.*blackman(length(x));
        
        % Compute power spectrum and normalize to probability density
        % function
        N = length(x);
        xdft = fft(x);
        xdft = xdft(1:(N/2+1));
        psdx = abs(xdft).^2;
        pdfx= psdx./sum(psdx) ; % Normalize to get a probability density function
        
        % Entropy Calculation
        for j = 1:length(freq)-freq_win_sz
           aux1=pdfx(j:j+freq_win_sz);
           entropy_eeg2pain(k,j)=-sum(aux1.*log2(aux1 + 1e-12))/log2(length(aux1));
        end
        
     end     
     
     for k=1:temp_win_step_sz:length(EEG2NoPain)-temp_win_sz
        x=EEG2NoPain(k:k+temp_win_sz);
        x=x.*blackman(length(x));
        
        % Compute power spectrum and normalize to probability density
        % function
        N = length(x);
        xdft = fft(x);
        xdft = xdft(1:(N/2+1));
        psdx = abs(xdft).^2;
        pdfx= psdx./sum(psdx) ; % Normalize to get a probability density function
        
        % Entropy Calculation
        for j = 1:length(freq)-freq_win_sz
           aux1=pdfx(j:j+freq_win_sz);
           entropy_eeg2nopain(k,j)=-sum(aux1.*log2(aux1 + 1e-12))/log2(length(aux1));
        end
        
        
     end          
     
 ppgpain(u,:)=mean(entropy_ppgpain);
 ppgnopain(u,:)=mean(entropy_ppgnopain); 
 eeg1pain(u,:)=mean(entropy_eeg1pain);
 eeg1nopain(u,:)=mean(entropy_eeg1nopain);
 eeg2pain(u,:)=mean(entropy_eeg2pain);
 eeg2nopain(u,:)=mean(entropy_eeg2nopain);
 
 end
 
 %% Wilcoxon signed rank test

% p-values calc, frec-wise and ln(p) vs frec plot
for k=1:size(ppgpain,2)
    [pPPG,h,stats]=signrank(ppgpain(:,k), ppgnopain(:,k));
    PPPG(k) = pPPG;
end

for k=1:size(eeg1pain,2)
    [pEEG1,h,stats]=signrank(eeg1pain(:,k), eeg1nopain(:,k));
    PEEG1(k) = pEEG1;
end

for k=1:size(eeg2pain,2)
    [pEEG2,h,stats]=signrank(eeg2pain(:,k), eeg2nopain(:,k));
    PEEG2(k) = pEEG2;
end

%% Plots
figure;

y=1:245;
a = log10(PPPG);
plot(y, a);
xlabel('Frequency onset (Hz)');
ylabel('ln(p-val)');
%title('PPG Wilcoxon signed rank test Pain VS No Pain');

hold all;
b = log10(PEEG1);
plot(y, b);
xlabel('Frequency (Hz)');
ylabel('Discriminability between conditions C1 and C2 (log10 pval)');
%title('EEG1 Wilcoxon signed rank test Pain VS No Pain');

hold all;
c = log10(PEEG2);
plot(y, c);
%title('EEG2 Wilcoxon signed rank test Pain VS No Pain');


legend('PPG','EEG1','EEG2');
set(legend);
%title('GSR Wilcoxon signed rank test Pain VS No Pain');

%toc

