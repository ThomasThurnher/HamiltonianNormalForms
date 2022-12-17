function [omega] = omega_bb(I,gamma,omega0)

omega = imag(omega0)*ones(size(I));

for j = 1:length(gamma)
    omega = omega + imag(gamma(j)) * (I.^j);
end

end
