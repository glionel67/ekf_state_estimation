function q = quaternionProduct(q1, q2)
%quaternionProduct: product of two quaternions
% !WARNING! not commutative: q1 x q2 != q2 x q1
% Args:
% - q1 (vector [4x1]): 1st quaternion to multiply
% - q2 (vector [4x1]): 2nd quaternion to multiply
% Return:
% - q (vector [4x1]): quaternion product

assert(length(q1) == 4, 'quaternionNormalize: q1 must be a [4x1] vector!');
assert(length(q2) == 4, 'quaternionNormalize: q2 must be a [4x1] vector!');

qw = q1(1)*q2(1) - q1(2)*q2(2) - q1(3)*q2(3) - q1(4)*q2(4);
qx = q1(1)*q2(2) + q1(2)*q2(1) + q1(3)*q2(4) - q1(4)*q2(3);
qy = q1(1)*q2(3) - q1(2)*q2(4) + q1(3)*q2(1) + q1(4)*q2(2);
qz = q1(1)*q2(4) + q1(2)*q2(3) - q1(3)*q2(2) + q1(4)*q2(1);

q = [qw; qx; qy; qz];
end
