function pu = normal_distribution(x, mu, Sigma)
    pu = (1/(2*pi*Sigma))*exp(-((x - mu)^2)/(2*Sigma));
end