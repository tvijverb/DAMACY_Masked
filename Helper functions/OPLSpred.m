function [tnew,Tnew_o,yhat, Enew] = OPLSpred(Xnew,P_o,W_o,w,q,LV)

% Johan Westerhuis
% Biosystems Data Analysis
% University of AMsterdam

% JAW: 20 juli 2007


Enew = Xnew;
Tnew_o = [];

for lv = 1:LV
    w_o = W_o(:,lv); 
    p_o = P_o(:,lv);
    tnew_o = w_o'*w_o\Enew*w_o;
    Tnew_o = [Tnew_o tnew_o];
    Enew = Enew - tnew_o*p_o';
end

tnew = w'*w\Enew*w;
yhat = tnew*q;