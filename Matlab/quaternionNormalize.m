function q = quaternionNormalize(q)
%quaternionNormalize: normalize a quaternion q
% Args:
% - q (vector [4x1]): quaternion to normalize
% Return:
% - q (vector [4x1]): normalized quaternion

assert(length(q) == 4, 'quaternionNormalize: q must be a [4x1] vector!');

%qNorm = norm(q, 2);
qNorm = sqrt(q(1)^2 + q(2)^2 + q(3)^2 + q(4)^2);

if (qNorm < 1e-9) % Avoid zero division!
    q = [1.0; 0.0; 0.0; 0.0];
else
    q = q ./ qNorm;
end

end
