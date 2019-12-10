clear all
clc

%% Unispectral analysis

% Load data
filePattern=fullfile('/path/to/data/*.txt');
files=dir(filePattern);

% Initialize variables
ppgpain=zeros(36,245);
ppgnopain=zeros(36,245);
eeg1pain=zeros(36,245);
eeg1nopain=zeros(36,245);
eeg2pain=zeros(36,245);
eeg2nopain=zeros(36,245);


 for u=1:length(files) % for each patient
     filename=files(u);
     name=filename.name;
     X=load(name);
     PPG=X(:,1);
     EEG1=X(:,3);
     EEG2=X(:,4);
     % GSR=X(:,2);
     PULS=X(:,5);
     win=500; % FFT calc window size
     win2=5; % Entropies calc window size
     
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
        
    % Entropies calc
    
    for k=1:length(PPGPain)-win
        x=PPGPain(k:k+win);
        f=fft(x,length(x));
        f=f(round(1:end/2));
        pow=abs(f).^2;
        pf1=(pow)./sum(pow);
        for j = 1:250-win2
            aux1=pf1(j:j+win2);
            entropy_ppgpain(k,j)=sum(aux1.*log(1./aux1))./log(length(aux1));
        end
    end
    
      
    for k=1:length(PPGNoPain)-win
        x=PPGNoPain(k:k+win);
        f=fft(x,length(x));
        f=f(round(1:end/2));
        pow=abs(f).^2;
        pf1=(pow)./sum(pow);
        for j = 1:250-win2
            aux1=pf1(j:j+win2);
            entropy_ppgnopain(k,j)=sum(aux1.*log(1./aux1))./log(length(aux1));
        end
    end
    
    
     for k=1:length(EEG1Pain)-win
        x=EEG1Pain(k:k+win);
        f=fft(x,length(x));
        f=f(round(1:end/2));
        pow=abs(f).^2;
        pf1=(pow)./sum(pow);
        for j = 1:250-win2
            aux1=pf1(j:j+win2);
            entropy_eeg1pain(k,j)=sum(aux1.*log(1./aux1))./log(length(aux1));
        end
     end
    
     for k=1:length(EEG1NoPain)-win
        x=EEG1NoPain(k:k+win);
        f=fft(x,length(x));
        f=f(round(1:end/2));
        pow=abs(f).^2;
        pf1=(pow)./sum(pow);
        for j = 1:250-win2
            aux1=pf1(j:j+win2);
            entropy_eeg1nopain(k,j)=sum(aux1.*log(1./aux1))./log(length(aux1));
        end
     end
    
     for k=1:length(EEG2Pain)-win
        x=EEG2Pain(k:k+win);
        f=fft(x,length(x));
        f=f(round(1:end/2));
        pow=abs(f).^2;
        pf1=(pow)./sum(pow);
        for j = 1:250-win2
            aux1=pf1(j:j+win2);
            entropy_eeg2pain(k,j)=sum(aux1.*log(1./aux1))./log(length(aux1));
        end
     end     
     
     for k=1:length(EEG2NoPain)-win
        x=EEG2NoPain(k:k+win);
        f=fft(x,length(x));
        f=f(round(1:end/2));
        pow=abs(f).^2;
        pf1=(pow)./sum(pow);
        for j = 1:250-win2
            aux1=pf1(j:j+win2);
            entropy_eeg2nopain(k,j)=sum(aux1.*log(1./aux1))./log(length(aux1));
        end
     end          
     
 ppgpain(u,:)=mean(entropy_ppgpain);
 ppgnopain(u,:)=mean(entropy_ppgnopain); 
 eeg1pain(u,:)=mean(entropy_eeg1pain);
 eeg1nopain(u,:)=mean(entropy_eeg1nopain);
 eeg2pain(u,:)=mean(entropy_eeg2pain);
 eeg2nopain(u,:)=mean(entropy_eeg2nopain);
 
 end

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

figure;
y=(1:245);
plot(y, log10(PPPG));
xlabel('Frequency onset (Hz)');
ylabel('ln(p-val)');
%title('PPG Wilcoxon signed rank test Pain VS No Pain');

hold all;
y=(1:245);
plot(y, log10(PEEG1));
xlabel('Freq');
ylabel('p-val');
%title('EEG1 Wilcoxon signed rank test Pain VS No Pain');

hold all;
y=(1:245);
plot(y, log10(PEEG2));
xlabel('Freq');
ylabel('p-val');
%title('EEG2 Wilcoxon signed rank test Pain VS No Pain');


legend('PPG','EEG1','EEG2');
set(legend);
%title('GSR Wilcoxon signed rank test Pain VS No Pain');
