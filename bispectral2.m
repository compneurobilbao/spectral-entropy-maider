% clear all
clc

tic

% Unispectral analysis
Fs = 500;  % Sampling frequency  
temp_win_sz=500; % Temporal window size (temp_win_sz/Fs = seconds of )
temp_win_step_sz=1; % Number of moving samples between one temporal window and the next one

freq = 0:Fs/temp_win_sz:Fs/2;

freq_win_sz=5; % Entropies calc window size
  

% Load data
filePattern=fullfile('/path/to/data/*.txt');
files=dir(filePattern);
num_users=length(files);

% Initialize variables

jointentropyPainPE1=zeros(num_users,length(freq)-freq_win_sz,length(freq)-freq_win_sz);
jointentropyNoPainPE1=zeros(num_users,length(freq)-freq_win_sz,length(freq)-freq_win_sz);
jointentropyPainPE2=zeros(num_users,length(freq)-freq_win_sz,length(freq)-freq_win_sz);
jointentropyNoPainPE2=zeros(num_users,length(freq)-freq_win_sz,length(freq)-freq_win_sz);
jointentropyPainE1E2=zeros(num_users,length(freq)-freq_win_sz,length(freq)-freq_win_sz);
jointentropyNoPainE1E2=zeros(num_users,length(freq)-freq_win_sz,length(freq)-freq_win_sz);




% for each user
 for u=1:num_users
          
     filename=files(u);
     name=filename.name;
     X=load(name);
     PPG=X(:,1);
     EEG1=X(:,3);
     EEG2=X(:,4);
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
    PPGPain=PPG(n_startpain:m,:);
    PPGNoPain=PPG(m_startpain:end,:);
    EEG1Pain=EEG1(n_startpain:m,:);
    EEG1NoPain=EEG1(m_startpain:end,:);
    EEG2Pain=EEG2(n_startpain:m,:);
    EEG2NoPain=EEG2(m_startpain:end,:);
        
    %% Entropies calc PPGEEG1
    N = temp_win_sz;
   
    for k=1:temp_win_step_sz:length(PPGPain)-temp_win_sz0 
        data_pain=[PPGPain(k:k+temp_win_sz).*blackman(length(temp_win_sz)) EEG1Pain(k:k+temp_win_sz).*blackman(length(temp_win_sz))];
        data_nopain=[PPGNoPain(k:k+temp_win_sz).*blackman(length(temp_win_sz)) EEG1NoPain(k:k+temp_win_sz).*blackman(length(temp_win_sz))];

        % Compute power spectrum and normalize to probability density
        % function
        xdft = fft2(data_pain, length(data_pain), 2);
        xdft = xdft(1:N/2+1,:);
        psdx = abs(xdft).^2;
        pdfx_pain= psdx./sum(psdx(:)) ; % Normalize to get a probability density function

        xdft = fft2(data_nopain, length(data_nopain), 2);
        xdft = xdft(1:N/2+1,:);
        psdx = abs(xdft).^2;
        pdfx_nopain= psdx./sum(psdx(:)) ; % Normalize to get a probability density function


        
        % Entropy Calculation  
        for k1=1:length(freq)-freq_win_sz
            for k2=1:length(freq)-freq_win_sz
                 aux1=pdfx_pain(k1:k1+freq_win_sz,1);
                 aux2=pdfx_pain(k2:k2+freq_win_sz,2);
                 jointentropyEEGPPGPain(k,k1,k2)=-sum(aux1.*log2(aux1 + 1e-12))/log2(length(aux1))-sum(aux2.*log2(aux2 + 1e-12))/log2(length(aux2));
                 
                 aux1=pdfx_nopain(k1:k1+freq_win_sz,1);
                 aux2=pdfx_nopain(k2:k2+freq_win_sz,2);
                jointentropyEEGPPGNoPain(k,k1,k2)=-sum(aux1.*log2(aux1 + 1e-12))/log2(length(aux1))-sum(aux2.*log2(aux2 + 1e-12))/log2(length(aux2));                
            end
        end
    end

   
   jointentropyPainPE1(u,:,:)=squeeze(mean(jointentropyEEGPPGPain));
   jointentropyNoPainPE1(u,:,:)=squeeze(mean(jointentropyEEGPPGNoPain));
   clear jointentropyEEG1PPGPain jointentropyEEG1PPGNoPain;
   

   
   %% PPGEEG2
   
   for k=1:temp_win_step_sz:length(PPGPain)-temp_win_sz
       
        data_pain=[PPGPain(k:k+temp_win_sz).*blackman(length(temp_win_sz)) EEG2Pain(k:k+temp_win_sz).*blackman(length(temp_win_sz))];
        data_nopain=[PPGNoPain(k:k+temp_win_sz).*blackman(length(temp_win_sz)) EEG2NoPain(k:k+temp_win_sz).*blackman(length(temp_win_sz))];

        % Compute power spectrum and normalize to probability density
        % function
        xdft = fft2(data_pain, length(data_pain), 2);
        xdft = xdft(1:N/2+1,:);
        psdx = abs(xdft).^2;
        pdfx_pain= psdx./sum(psdx(:)) ; % Normalize to get a probability density function

        xdft = fft2(data_nopain, length(data_nopain), 2);
        xdft = xdft(1:N/2+1,:);
        psdx = abs(xdft).^2;
        pdfx_nopain= psdx./sum(psdx(:)) ; % Normalize to get a probability density function

        % Entropy Calculation  
        for k1=1:length(freq)-freq_win_sz
            for k2=1:length(freq)-freq_win_sz
                 aux1=pdfx_pain(k1:k1+freq_win_sz,1);
                 aux2=pdfx_pain(k2:k2+freq_win_sz,2);
                 jointentropyEEGPPGPain(k,k1,k2)=-sum(aux1.*log2(aux1 + 1e-12))/log2(length(aux1))-sum(aux2.*log2(aux2 + 1e-12))/log2(length(aux2));
                
                 aux1=pdfx_nopain(k1:k1+freq_win_sz,1);
                 aux2=pdfx_nopain(k2:k2+freq_win_sz,2);
                 jointentropyEEGPPGNoPain(k,k1,k2)=-sum(aux1.*log2(aux1 + 1e-12))/log2(length(aux1))-sum(aux2.*log2(aux2 + 1e-12))/log2(length(aux2));                
           end
        end
   end


   jointentropyPainPE2(u,:,:)=squeeze(mean(jointentropyEEGPPGPain));
   jointentropyNoPainPE2(u,:,:)=squeeze(mean(jointentropyEEGPPGNoPain));
   clear jointentropyEEG2PPGPain jointentropyEEG2PPGNoPain;
   
   %% EGG1 EGG2
   
   for k=1:temp_win_step_sz:length(PPGPain)-temp_win_sz
       
        data_pain=[EEG1Pain(k:k+temp_win_sz).*blackman(length(temp_win_sz)) EEG2Pain(k:k+temp_win_sz).*blackman(length(temp_win_sz))];
        data_nopain=[EEG1NoPain(k:k+temp_win_sz).*blackman(length(temp_win_sz)) EEG2NoPain(k:k+temp_win_sz).*blackman(length(temp_win_sz))];

        % Compute power spectrum and normalize to probability density
        % function
        xdft = fft2(data_pain, length(data_pain), 2);
        xdft = xdft(1:N/2+1,:);
        psdx = abs(xdft).^2;
        pdfx_pain= psdx./sum(psdx(:)) ; % Normalize to get a probability density function

        xdft = fft2(data_nopain, length(data_nopain), 2);
        xdft = xdft(1:N/2+1,:);
        psdx = abs(xdft).^2;
        pdfx_nopain= psdx./sum(psdx(:)) ; % Normalize to get a probability density function

        % Entropy Calculation  
        for k1=1:length(freq)-freq_win_sz
            for k2=1:length(freq)-freq_win_sz
                 aux1=pdfx_pain(k1:k1+freq_win_sz,1);
                 aux2=pdfx_pain(k2:k2+freq_win_sz,2);
                 jointentropyEEGPain(k,k1,k2)=-sum(aux1.*log2(aux1 + 1e-12))/log2(length(aux1))-sum(aux2.*log2(aux2 + 1e-12))/log2(length(aux2));
                
                 aux1=pdfx_nopain(k1:k1+freq_win_sz,1);
                 aux2=pdfx_nopain(k2:k2+freq_win_sz,2);
                 jointentropyEEGNoPain(k,k1,k2)=-sum(aux1.*log2(aux1 + 1e-12))/log2(length(aux1))-sum(aux2.*log2(aux2 + 1e-12))/log2(length(aux2));
                 
             end
        end
  
   end

   jointentropyPainE1E2(u,:,:)=squeeze(mean(jointentropyEEGPain));
   jointentropyNoPainE1E2(u,:,:)=squeeze(mean(jointentropyEEGNoPain));
   clear jointentropyEEGsPain jointentropyEEGsNoPain;

 end
 
 

%% p-values calc, frec-wise and ln(p) vs frec plot

for k1=1:length(freq)-freq_win_sz
    for k2=1:length(freq)-freq_win_sz
        p1(k1,k2)=signrank(jointentropyPainPE1(:,k1,k2), jointentropyNoPainPE1(:,k1,k2));
        p2(k1,k2)=signrank(jointentropyPainPE2(:,k1,k2), jointentropyNoPainPE2(:,k1,k2));
        p3(k1,k2)=signrank(jointentropyPainE1E2(:,k1,k2), jointentropyNoPainE1E2(:,k1,k2));
        
    end
end

% plots
subplot(1,3,1)
n=log10(p1);
imagesc(n);
title('PPG + EEG1');
colorbar
caxis([-5 0]);
axis square
ylabel("freq. of PPG (Hz)");
xlabel("freq. of EEG1 (Hz)");

subplot(1,3,2)
m=log10(p2);
imagesc(m);
title('PPG + EEG2');
colorbar
caxis([-5 0]);
axis square
ylabel("freq. of PPG (Hz)");
xlabel("freq. of EEG2 (Hz)");


subplot(1,3,3)
h=log10(p3);
imagesc(h);
title('EEG1 + EEG2');
colorbar
caxis([-5 0]);
axis square
ylabel("freq. of EEG1 (Hz)");
xlabel("freq. of EEG2 (Hz)");

toc
