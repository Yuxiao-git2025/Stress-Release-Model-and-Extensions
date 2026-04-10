function [J,A,B,C]=LSRM_UnpackingParam(params)
J = round(sqrt(numel(params)+1)-1);
% J=round((-2+sqrt(4+4*n))/2);
A = params(1:J); A=A(:);               % Jx1
B = params(J+1:2*J); B=B(:);           % Jx1
C = reshape(params(2*J+1:end), [J, J]);
end