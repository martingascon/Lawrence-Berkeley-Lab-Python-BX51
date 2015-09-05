
# Gaussian function
def gauss_func(x, *p):
    A, mu, sigma = p
    return (A/(sqrt(2*pi)*sigma)) * exp(-(x-mu)**2 /(2*sigma**2))
