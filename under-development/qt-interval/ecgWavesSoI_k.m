function [P, Q, R, S, T] = ecgWavesSoI_k(ecg, Rpeak, RRf)

% 'knot based' version.
% This function appriximately finds peaks and borders of the Q and T waves 
% in the input ecg.   


% Inputs:
% ecg: single channel ecg signal,
% R: R peak positions,
% RRf: R-peak referenced flag; if true, R peak is the reference of all the
% output points. if false, they are based on the samples index.

% Outputs:
% P, Q, R, S, T: Each of the outputs is a matrix, in which the first, second, and third 
% columns respectively contain onset, peak and offset of the P, Q, R, S, T waves. 


% Davood Fattahi, 5/5/2021
% fattahi.d@gmail.com




if length(Rpeak)==length(ecg)
    Rpeak = find(Rpeak);
end

ecg=ecg(:); Rpeak=Rpeak(:);
RR=Rpeak(2:end)-Rpeak(1:end-1); RR = [median(RR); RR; median(RR)];
P=nan(length(Rpeak),3); Q=nan(length(Rpeak),3);
R=nan(length(Rpeak),3); S=nan(length(Rpeak),3);
T=nan(length(Rpeak),3); 
% Uafp=nan(length(Rpeak),3);


for i= 1: length(RR)-1
    %% Q
    s=Rpeak(i)-floor(.4*RR(i));
    e=Rpeak(i);
    if s<1
        s=1;
    end
    if e>length(ecg)
        e=length(ecg);
    end
    k = findknots(ecg(s:e),3,'Value') + s;
    Q(i,2)=k(2);
    Q(i,1)=Q(i,2)-ceil(.050*RR(i)); % it is also the iso-point
    [~,I]=min(abs(ecg(Q(i,2):Rpeak(i))-ecg(Q(i,1)))); % the nearest point to iso-point
    Q(i,3)=Q(i,2)+I-1;
    
    %% P
    s=Rpeak(i)-floor(.4*RR(i));
    e=Q(i,1);
    if s<1
        s=1;
    end
    if e>length(ecg)
        e=length(ecg);
    end
    k = findknots(ecg(s:e),3,'Value') + s;  
    P(i,:)=k;
    
    %% S
    s=Rpeak(i);
    e=Rpeak(i)+ceil(.15*RR(i+1));
    if s<1
        s=1;
    end
    if e>length(ecg)
        e=length(ecg);
    end
    k = findknots(ecg(s:e),3,'Value') + s;
    S(i,2)=k(2);
    S(i,3)=S(i,2)+ceil(.050*RR(i));
    [~,I]=min(abs(ecg(s:S(i,2))-ecg(S(i,3)))); % the nearest point to ST level
    S(i,1)=s+I-1;
    
    %% R
    R(i,1)=Q(i,3); R(i,2)=Rpeak(i); R(i,3) = S(i,1);
    %% T
    s=Rpeak(i)+floor(.1*RR(i+1));
    e=Rpeak(i)+floor(.6*RR(i+1));
    if s<1
        s=1;
    end
    if e>length(ecg)
        e=length(ecg);
    end
    k = findknots(ecg(s:e),3,'Value') + s;
    T(i,2)=k(2);
    T(i,3)=T(i,2)+floor(.25*RR(i));
    T(i,1)=T(i,2)-floor(.20*RR(i));
    
    %% 
end


P(P<1|P>length(ecg))=nan;
Q(Q<1|Q>length(ecg))=nan;
R(R<1|R>length(ecg))=nan;
S(S<1|S>length(ecg))=nan;
T(T<1|T>length(ecg))=nan;


if RRf
    P=P-Rpeak; Q=Q-Rpeak; R=R-Rpeak; S=S-Rpeak; T=T-Rpeak;
end

end

