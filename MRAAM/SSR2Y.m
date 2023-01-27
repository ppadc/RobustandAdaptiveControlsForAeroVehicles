function [Acl, Bcl, Ccl, Dcl, Bcl2, Dcl2] = SSR2Y(Ap, Bp, Cp, Dp, Ac, ...
    Bc1, Bc2, Cc, Dc1, Dc2)

Iu = eye(size(Dc1,1));

Z = Iu - Dc1*Dp;
Zinv = Z\Iu;

Acl = [Ap + Bp*Zinv*Dc1*Cp, Bp*Zinv*Cc; Bc1*Cp + Bc1*Dp*Zinv*Dc1*Cp, ...
    Ac + Bc1*Dp*Zinv*Cc];

Bcl = [ Bp*Zinv*Dc2; Bc1*Dp*Zinv*Dc2 + Bc2];

Ccl = [Cp + Dp*Zinv*Dc1*Cp, Dp*Zinv*Cc];

Dcl = Dp*Zinv*Dc2;

% These matrices are specific to an additive noise being placed on the
% sensor output y = Cx + Du + v, where v is the noise term.
Bcl2 = [ Bp*Zinv*Dc1; Bc1*Dp*Zinv*Dc1 + Bc1];
Dcl2 = Dp*Zinv*Dc1 + eye(size(Dp,1));