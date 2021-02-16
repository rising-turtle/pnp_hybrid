function [error] = EvaluateSampsonError(point1,point2,F)
nbrpoints = size(point1,2);

if size(point1,1) == 2
    point1 = [point1;ones(1,nbrpoints)];
    point2 = [point2;ones(1,nbrpoints)];
end

Mat1 = F*point2;
Mat2 = F'*point1;
deno = Mat1(1,:).^2 + Mat1(2,:).^2 + Mat2(1,:).^2 + Mat2(2,:).^2;

F = F.';
vecF = F(:);
u1 = point1(1,:).'; v1 = point1(2,:).';
u2 = point2(1,:).'; v2 = point2(2,:).';

M = [u1.*u2 u1.*v2 u1 v1.*u2 v1.*v2 v1 u2 v2 ones(nbrpoints,1)];
nomi = (M*vecF).^2;

error = sum((nomi.')./deno);
%error = sqrt(error/nbrpoints);
end