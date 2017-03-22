function [w,q,P_o,W_o] = OPLS(X,y,LV)

% Johan Westerhuis
% Biosystems Data Analysis
% University of AMsterdam

% LV orthogonal components are calculated after which 1 predictive
% component is calculated

% JAW: 20 juli 2007

T_o = zeros(length(X(:,1)),LV);
P_o = zeros(length(X(1,:)),LV);
W_o = zeros(length(X(1,:)),LV);

E = X;
w = ((y'*y)\y'*E)';  
w = w  / sqrt(w'*w);

for lv = 1:LV;    
    t = (w'*w)\E*w;
    % c = (inv(t'*t)*t'*y)'; 
    % u = (inv*c'*c)*y*c;
    % These are two strange lines. For single y c is a scaler (so why
    % transposing it). Furthermore, c and u are not used anymore
    p = ((t'*t)\t'*E)';

    w_o = p - ((w'*w)\w'*p)*w;
    w_o = w_o / sqrt(w_o'*w_o);
    t_o = (w_o'*w_o)\E*w_o;
    p_o = ((t_o'*t_o)\t_o'*E)';

    E = E - t_o*p_o';

    T_o(:,lv) = t_o;
    P_o(:,lv) = p_o;
    W_o(:,lv) = w_o;
end

w = ((y'*y)\y'*E)';
w = w / sqrt(w'*w);
t = (w'*w)\E*w;
q = ((t'*t)\t'*y)';
