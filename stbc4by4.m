%%%%%%%%%%%%% 4 X 4 STBC %%%%%%%%%%%%%%%
 
clc
clear all;
close all;
%parameter stup
N = 4*10^4; %number of bit to process
M=4; %size of signal constellation
x = randi([0 3],1,N); %generating random number

% signal modulation
xmod=qammod(x,4);
xmod=reshape(xmod,4,N/4);
xmod=kron(xmod,[1,1,1,1]);%Matrix formed by taking all possible products
%between the elements of xmod and the matrix of 1's
 
%generating 4x4 matrix and zero forcing
for i=1:16
    h(i,:)=1/sqrt(2)*(randn(1,N/4) + 1i*randn(1,N/4)); 
end
H=reshape(h,4,N);
y=reshape(sum(H.*xmod,1),4,N/4);
H=reshape(h,4,4,N/4);
snr=linspace(0,6,7);
ser=zeros(1,length(snr));
for ii=1:length(snr)
    N1=1/sqrt(2)*(randn(1,N)+1i*randn(1,N));
    N1=reshape(N1,4,N/4); 
    ynoisy=y+10^(-(snr(ii)-10*log10(16))/20)*N1;
    ynoisy=awgn(y,snr(ii),'measured'); %adding white gaussian noise
    ynoisy=reshape(ynoisy,4,1,N/4);
    B=[];recvd=[];
    for kk=1:N/4
        Heq=transpose(H(:,:,kk));
        B=pinv(Heq);
        recvd=[recvd,B*ynoisy(:,:,kk)];
    end
    finy=qamdemod(reshape(recvd,1,N),4); %signal demodulation
    [num ty]=symerr(x,finy);
    ser(ii)=ty;
end
semilogy(snr,ser,'r-*'); %result presentation
grid on;hold on;
title('Plot of Symbol error rate for 4X4(4-QAM) System','FontSize',12);
legend('sim (nTx=4, nRx=4, Uncoded(4-QAM))','location','southwest');
xlabel('SNR(dB) ---->','Color','k','FontSize',11);Ylabel('Symbol Error rate ---->','Color','k','FontSize',11);
