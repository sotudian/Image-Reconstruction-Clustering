function [V, U, T, E,Objective_Function_PFCM,Objective_Function_FCM,Objective_Function_PCM] = PFCM_clustering (X, c,m,eta,Cf,Cp,w)
w = w(:);
n = size(X, 1);
p = size(X, 2);

max_iter = 100;		% Max. iteration
term_thr = 1e-4;		% Termination threshold
display = 1;		% Display info or not
E = zeros(max_iter, 1);	% Array for termination measure values
Objective_Function_PCM=zeros(max_iter, 1);
Objective_Function_PFCM=zeros(max_iter, 1);
Objective_Function_FCM=zeros(max_iter, 1);
V = rand(c, p);

U = zeros (c, n);
T = zeros (c, n);

% Main loop
for i = 1:max_iter,
    
    
% fill the distance matrix
dist = Distance_Function (V, X); 
% calculate new U, suppose m != 1
tmp = dist.^(-2/(m-1));      
U = tmp./(ones(c, 1)*sum(tmp));

% Correct the situation of "singularity" (one of the data points is
% exactly the same as one of the cluster centers).
si = find (tmp == Inf);
U(si) = 1;
if (size (si, 1) ~= 0)
    display ('FPCMC, Warning: Singularity occured and corrected.');
end


% Claculate new T

tmp = (Cp.*(dist .^ 2)) ./ ( w * ones (1, n));
tmp = tmp.^(1/(eta-1));
T = 1 ./ (1 + tmp);

% Correct the situation of "singularity" (one of the data points is
% exactly the same as one of the cluster centers).
T(si) = 1; % Do more later


% % Check constraint
% tmp = find ((sum (T') - ones (1, c)) > 0.0001);
% if (size(tmp,2) ~= 0)
%     display ('FPCMC, Warning: Constraint for T is not hold.');
% end

% Scale T
% maxU =  max(max(U));
% maxT =  max(max(T));
% if maxT < .5
%    ScF = .5 / maxT;
%    T = T * ScF;
% end

%   objective function FCM
PU1 = U.^m;
Objective_Function_FCM(i) = sum(sum((dist.^2).*PU1)); 
%   objective function PCM
PU2 = T.^eta;
Objective_Function_PCM(i) = sum(sum((dist.^2).*PU2))+sum(w.*sum(((1-T).^eta),2)); 
%   objective function PFCM
Objective_Function_PFCM(i) = Cf*Objective_Function_FCM(i)  +  Cp*sum(sum((dist.^2).*PU2))  +  sum(w.*sum(((1-T).^eta),2));


V_old = V;
% new center
Us = Cf.*(U.^m);
Ts = Cp.*(T.^eta);
V = ((Us+Ts)*X) ./ ((ones(p, 1)*sum((Us+Ts)'))'); 

E(i) = norm (V - V_old, 1);


  
    if display, 
		fprintf('Iteration count = %d, Termination measure value = %f\n', i, E(i));
	end

    % check termination condition
	if E(i) <= term_thr, break; end,
end

% iter_n = i;	% Actual number of iterations 
% E(iter_n+1:max_iter) = [];
