%%%%%%%%%%%%% 128 X 128 STBC %%%%%%%%%%%%%%%
 
clc
clear all;
close all;
%parameter stup
N = 128*10;
M=4;
x = randi([0 3],1,N); %generating random number
%modulation
xmod=qammod(x,4);
xmod=reshape(xmod,128,N/128);
xmod=kron(xmod,[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1]);
 
%generating 4x4 matrix
for i=1:(128*128)
    h(i,:)=1/sqrt(2)*(randn(1,N/128) + 1i*randn(1,N/128)); 
end
H=reshape(h,128,N);
y=reshape(sum(H.*xmod,1),128,N/128);
H=reshape(h,128,128,N/128);
snr=linspace(0,6,7);
ser=zeros(1,length(snr));
for ii=1:length(snr)
    N1=1/sqrt(2)*(randn(1,N)+1i*randn(1,N));
    N1=reshape(N1,128,N/128);
    ynoisy=y+10^(-(snr(ii)-10*log10(128*128))/20)*N1;
    ynoisy=awgn(y,snr(ii),'measured'); %adding white gaussian noise
    ynoisy=reshape(ynoisy,128,1,N/128);
    B=[];recvd=[];
    for kk=1:N/128
        Heq=transpose(H(:,:,kk));
        B=pinv(Heq);
        recvd=[recvd,B*ynoisy(:,:,kk)];
    end
    finy=qamdemod(reshape(recvd,1,N),128); %demodulation
    [num ty]=symerr(x,finy);
    ser(ii)=ty;
end
semilogy(snr,ser,'r-*');
grid on;hold on;
title('Plot of Symbol error rate for 128X128(4-QAM) System','FontSize',12);
legend('sim (nTx=128, nRx=128, Uncoded(4-QAM))','location','southwest');
xlabel('SNR(dB) ---->','Color','k','FontSize',11);Ylabel('Symbol Error rate ---->','Color','k','FontSize',11);
