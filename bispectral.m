clear all
clc

%% Bispectral analysis

filePattern=fullfile('/path/to/data/*.txt');
files=dir(filePattern);

for u=1:length(files) % for each patient

     filename=files(u);
     name=filename.name;
     X=load(name);
     PPG=X(:,1);
     EEG1=X(:,3);
     EEG2=X(:,4);
     PULS=X(:,5);
     win=500; % FFT calc window size
     win2=5; % Entropies calc window size
     
    % Pain / No Pain area limits
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
    

    % 55025 is the min amount of points for a subject
    
    PPGPain=PPGPaintot(1:55025,:);
    PPGNoPain=PPGNoPaintot(1:55025,:);
    EEG1Pain=EEG1Paintot(1:55025,:);
    EEG1NoPain=EEG1NoPaintot(1:55025,:);
    EEG2Pain=EEG2Paintot(1:55025,:);
    EEG2NoPain=EEG2NoPaintot(1:55025,:);


    data=[PPGPain EEG1Pain];
    f=fft2(data, length(data), 2);
    f=f(round(1:end/2),:);
    painpowPE1(:,:,u)=abs(f).^2;
    
    data=[PPGPain EEG2Pain];
    f=fft2(data, length(data), 2);
    f=f(round(1:end/2),:);
    painpowPE2(:,:,u)=abs(f).^2;
    
        
    data=[PPGNoPain EEG1NoPain];
    f=fft2(data, length(data), 2);
    f=f(round(1:end/2),:);
    nopainpowPE1(:,:,u)=abs(f).^2;
    
        
    data=[PPGNoPain EEG2NoPain];
    f=fft2(data, length(data), 2);
    f=f(round(1:end/2),:);
    nopainpowPE2(:,:,u)=abs(f).^2;
    
        
    data=[EEG1Pain EEG2Pain];
    f=fft2(data, length(data), 2);
    f=f(round(1:end/2),:);
    painpowE1E2(:,:,u)=abs(f).^2;
    
        
    data=[EEG1NoPain EEG2NoPain];
    f=fft2(data, length(data), 2);
    f=f(round(1:end/2),:);
    nopainpowE1E2(:,:,u)=abs(f).^2;
    
        
end

win=5;

jointentropyPE1=zeros(36,250-win,250-win);
jointentropyPE2=zeros(36,250-win,250-win);
jointentropyE1E2=zeros(36,250-win,250-win);

for k1=1:(250-win)
    for k2=1:(250-win)
        for i=1:36
            pf=squeeze(painpowPE1(:,:,i))./sum(sum(squeeze(painpowPE1(:,:,i))));

            aux1=pf(k1:k1+win,1);
            aux2=pf(k2:k2+win,2);

            jointentropyPE1(i,1,k1,k2)=sum(sum(-aux1.*log(aux1)))./log(numel(aux1))+sum(sum(-aux2.*log(aux2)))./log(numel(aux2));

            pf=squeeze(nopainpowPE1(:,:,i))./sum(sum(squeeze(nopainpowPE1(:,:,i))));

            aux1=pf(k1:k1+win,1);
            aux2=pf(k2:k2+win,2);

            jointentropyPE1(i,2,k1,k2)=sum(sum(-aux1.*log(aux1)))./log(numel(aux1))+sum(sum(-aux2.*log(aux2)))./log(numel(aux2));

        end
    end
end

for k1=1:(250-win)
    for k2=1:(250-win)
        for i=1:36
            pf=squeeze(painpowPE2(:,:,i))./sum(sum(squeeze(painpowPE2(:,:,i))));

            aux1=pf(k1:k1+win,1);
            aux2=pf(k2:k2+win,2);

            jointentropyPE2(i,1,k1,k2)=sum(sum(-aux1.*log(aux1)))./log(numel(aux1))+sum(sum(-aux2.*log(aux2)))./log(numel(aux2));

            pf=squeeze(nopainpowPE2(:,:,i))./sum(sum(squeeze(nopainpowPE2(:,:,i))));

            aux1=pf(k1:k1+win,1);
            aux2=pf(k2:k2+win,2);

            jointentropyPE2(i,2,k1,k2)=sum(sum(-aux1.*log(aux1)))./log(numel(aux1))+sum(sum(-aux2.*log(aux2)))./log(numel(aux2));

        end
    end
end

for k1=1:(250-win)
    for k2=1:(250-win)
        for i=1:36
            pf=squeeze(painpowE1E2(:,:,i))./sum(sum(squeeze(painpowE1E2(:,:,i))));

            aux1=pf(k1:k1+win,1);
            aux2=pf(k2:k2+win,2);

            jointentropyE1E2(i,1,k1,k2)=sum(sum(-aux1.*log(aux1)))./log(numel(aux1))+sum(sum(-aux2.*log(aux2)))./log(numel(aux2));

            pf=squeeze(nopainpowE1E2(:,:,i))./sum(sum(squeeze(nopainpowE1E2(:,:,i))));

            aux1=pf(k1:k1+win,1);
            aux2=pf(k2:k2+win,2);

            jointentropyE1E2(i,2,k1,k2)=sum(sum(-aux1.*log(aux1)))./log(numel(aux1))+sum(sum(-aux2.*log(aux2)))./log(numel(aux2));

        end
    end
end


p1=zeros(250-win, 250-win);
p2=zeros(250-win, 250-win);
p3=zeros(250-win, 250-win);



for k1=1:250-win
    for k2=1:250-win
        p1(k1,k2)=signrank(jointentropyPE1(1:36,1,k1,k2), jointentropyPE1(1:36,2,k1,k2));
        p2(k1,k2)=signrank(jointentropyPE2(1:36,1,k1,k2), jointentropyPE2(1:36,2,k1,k2));
        p3(k1,k2)=signrank(jointentropyE1E2(1:36,1,k1,k2), jointentropyE1E2(1:36,2,k1,k2));
        
    end
end



subplot(3,1,1)
n=log10(p1);
imagesc(n);
title('PPG + EEG1');
colorbar
caxis([-6 0]);

subplot(3,1,2)
m=log10(p2);
imagesc(m);
title('PPG + EEG2');
colorbar
caxis([-6 0]);


subplot(3,1,3)
h=log10(p3);
imagesc(h);
title('EEG1 + EEG2');
colorbar
caxis([-6 0]);

