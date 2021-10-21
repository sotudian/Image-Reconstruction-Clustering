function [V1,V2,UT1,UT2, E1,E2] = Interval_T2_PFCM_clustering (X, c,m1,m2,eta1,eta2,Cf,Cp,w)
w = w(:);
n = size(X, 1);
p = size(X, 2);
max_iter = 100;		% Max. iteration
term_thr = 1e-4;		% Termination threshold
display = 1;		% Display info or not

E1 = zeros(max_iter, 1);	% Array for termination measure values
E2 = zeros(max_iter, 1);	% Array for termination measure values
% Objective_Function_PFCM1=zeros(max_iter, 1);
% Objective_Function_PFCM2=zeros(max_iter, 1);
% Objective_Function_PFCM3=zeros(max_iter, 1);
% Objective_Function_PFCM4=zeros(max_iter, 1);


V1 = rand(c, p);
V2 = V1;

U1 = zeros (c, n);
T1 = zeros (c, n);
U2 = zeros (c, n);
T2 = zeros (c, n);

% Main loop
for i = 1:max_iter,
    
    
% fill the distance matrix
dist1 = Distance_Function (V1, X);
dist2 = Distance_Function (V2, X);

% calculate new U
tmp1 = dist1.^(-2/(m1-1));      
U_new1 = tmp1./(ones(c, 1)*sum(tmp1));
tmp2 = dist2.^(-2/(m2-1));      
U_new2 = tmp2./(ones(c, 1)*sum(tmp2));

U1=min(U_new1,U_new2);
U2=max(U_new1,U_new2);


% Correct the situation of "singularity" (one of the data points is
% exactly the same as one of the cluster centers).
si1 = find (tmp1 == Inf);
U1(si1) = 1;
si2 = find (tmp2 == Inf);
U2(si2) = 1;
if (size (si1, 1) ~= 0)
    disp('FPCMC, Warning: Singularity occured and corrected.');
end
if (size (si2, 1) ~= 0)
    disp('FPCMC, Warning: Singularity occured and corrected.');
end


% Claculate new T
tmp1 = (Cp.*(dist1 .^ 2)) ./ ( w * ones (1, n));
tmp1 = tmp1.^(1/(eta1-1));
T_new1= 1 ./ (1 + tmp1);
tmp2 = (Cp.*(dist2 .^ 2)) ./ ( w * ones (1, n));
tmp2 = tmp2.^(1/(eta2-1));
T_new2= 1 ./ (1 + tmp2);

T1=min(T_new1,T_new2);
T2=max(T_new1,T_new2);

% Correct the situation of "singularity" (one of the data points is
% exactly the same as one of the cluster centers).
T1(si1) = 1; % Do more later
T2(si2) = 1;




% % Check constraint
% tmp = find ((sum (T') - ones (1, c)) > 0.0001);
% if (size(tmp,2) ~= 0)
%     display ('FPCMC, Warning: Constraint for T is not hold.');
% end

% Scale T
%maxU =  max(max(U));
%maxT =  max(max(T));
%if maxT < .5
%    ScF = .5 / maxT;
%    T = T * ScF;
%end

% %   objective function PFCM
% PU1 = U1.^m1;
% PU1 = U1.^m1;
% PT1 = T1.^eta;
% Objective_Function_PFCM1(i) = Cf*sum(sum((dist.^2).*PU1))  +  Cp*sum(sum((dist.^2).*PT1))  +  sum(w.*sum(((1-T).^eta),2));
% Objective_Function_PFCM1(i) = Cf*sum(sum((dist.^2).*PU1))  +  Cp*sum(sum((dist.^2).*PT1))  +  sum(w.*sum(((1-T).^eta),2));
% Objective_Function_PFCM1(i) = Cf*sum(sum((dist.^2).*PU1))  +  Cp*sum(sum((dist.^2).*PT1))  +  sum(w.*sum(((1-T).^eta),2));
% Objective_Function_PFCM1(i) = Cf*sum(sum((dist.^2).*PU1))  +  Cp*sum(sum((dist.^2).*PT1))  +  sum(w.*sum(((1-T).^eta),2));

V_old1 = V1;
V_old2 = V2;
% Upper bound and Lower Bound Membership
U_T1=(Cf.*U1)+(Cp.*T1);
U_T2=(Cf.*U1)+(Cp.*T2);
U_T3=(Cf.*U2)+(Cp.*T1);
U_T4=(Cf.*U2)+(Cp.*T2);
% Lower Bound Membership
DUM1=min(U_T1,U_T2);
DUM2=min(DUM1,U_T3);
UT1=min(DUM2,U_T4);
% Upper bound Membership
DUM1=max(U_T1,U_T2);
DUM2=max(DUM1,U_T3);
UT2=max(DUM2,U_T4);
% % new center
% Us = Cf.*(U.^m);
% Ts = Cp.*(T.^eta);
% V = ((Us+Ts)*X) ./ ((ones(p, 1)*sum((Us+Ts)'))'); 


% Karni-Mendel Alg
for km=1:c
    % Hard Partitioning
UT_mean=(UT1+UT2)/2;
maxUT_mean = max(UT_mean);

    index_c = find(UT_mean(km, :) == maxUT_mean);
    Y= X(index_c,:);
       
    a=UT1(km,index_c);
    b=UT2(km,index_c);
    a=a(:);
    b=b(:);
    F=cat(2,a,b);
    if size(Y,1)<=1
    a1=Y;
    a2=km;
    a3=index_c;
    a4=UT_mean;
    a5=maxUT_mean;
    end
    
    [XLeft,XRight,L,R]=t2f_TR_KM(F,Y);
   
    V1(km,:)=XLeft';   
    V2(km,:)=XRight';
end





E1(i) = norm (V1 - V_old1, 1);
E2(i) = norm (V2 - V_old2, 1);

  
    if display, 
		fprintf('Iteration count = %d, Termination measure value => Lower= %f    Upper=%f\n', i, E1(i),E2(i));
	end

    % check termination condition
	if E1(i) <= term_thr, break; end,
    if E2(i) <= term_thr, break; end,
end


% iter_n = i;	% Actual number of iterations 
% E(iter_n+1:max_iter) = [];


end