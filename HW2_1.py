# Josh Popp
# jmp448
# CS 4775: Computational Genetics

"""
Script for computing h_{Y_n} given n and the probabilities of observing 1-6.
Arguments:
    n - integer representing the number of observations
    p - sequence of 6 float values summing to 1 representing pi_i, i in [1,6]

Outputs:
    h_<n>_dp.csv - a file the probabilities h_{Y_n}(y), y in [n,6n] in CSV
                   format, for use in later parts of question 1.
Example Usage:
    python 1c.py -n 50 -p .1 .2 .2 .2 .1 .2

    This would output a file called h_dp_50.csv.
"""

import argparse
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import math


"""
Computes h_{Y_n}(y) for all y in {n, ..., 6n} for a given n and pi.
Arguments:
    n - the number of random variables being summed
    pi - the probabilities of 1-6 for a single observation

Returns:
    a vector of probabilities for y in {n, ..., 6n}.
        index 0 should correspond to y = n, index 1 to n+1, etc.
"""


def h_Y(n, pi):
    # Initialize the n x 6n matrix
    h = np.zeros((n, 6*n))

    # Initial values
    # Probability of getting 1,...,6 in one roll is pi(1),...,pi(6)
    for i in range(6):
        h[0][i] = pi[i]

    for nj in range(1, n):
        for yi in range(1, 6 * n):
            if yi < 6:
                for k in range(yi):
                    h[nj][yi] += pi[k]*h[nj-1][yi-(k+1)]
            else:
                for k in range(6):
                    h[nj][yi] += pi[k]*h[nj-1][yi-(k+1)]
    return h[n-1][n-1:6*n]


""" Returns the minimum ten probabilities of a given array.
Arguments:
    n - the number of random variables being summed
    probs - the probabilities of [n, 6n]
Returns:
    a vector of the values of y in [n, 6n] that correspond to the minimum 10
        probabilities in probs priority is given to higher indices in case of
        ties.
"""


def min10(n, probs):
    # Create list of indices
    # Also keep track of minimum

    minima = list()
    minima.append(0)
    the_cut = 1
    for i in range(1, len(probs)):
        if probs[i] <= the_cut:
            j = 0
            while j < len(minima) and probs[i] > probs[minima[j]]:
                j += 1
            minima.insert(j, i)
            if len(minima) > 10:
                minima.pop()
                the_cut = probs[9]
    return minima


""" Returns the maximum ten probabilities of a given array.
Arguments:
    n - the number of random variables being summed
    probs - the probabilities of [n, 6n]
Returns:
    a vector of the values of y in [n, 6n] that correspond to the minimum 10
        probabilities in probs priority is given to higher indices in case of
        ties.
"""


def max10(n, probs):

    maxima = list()
    maxima.append(0)
    the_cut = 0
    for i in range(1,len(probs)):
        if probs[i] > the_cut:
            j = 0
            while j < len(maxima) and probs[i] < probs[maxima[j]]:
                j += 1
            maxima.insert(j, i)
            if len(maxima) > 10:
                maxima.pop(10)
                the_cut = probs[9]
    return maxima

""" Computes the mean and variance of an array of probabilities for [n, 6n].
Arguments:
    n - the value of n defining the range [n, 6n]
    probs - the array of probabilities
Returns:
    expect - the expectation
    variance - the variance
"""


def compute_stats(n, probs):
    # initialize:
    mean = 0.0
    variance = 0.0

    for i in range(len(probs)):
        mean += (i+n)*probs[i]

    for i in range(len(probs)):
        dif = i + n - mean
        variance += probs[i] * (dif ** 2)

    return mean, variance


def factorial(n):
    prod = 1
    while n > 1:
        prod *= n
        n -= 1
    return prod

def comb(n,k):
    prod = 1
    stop = n-k
    while n > stop:
        prod *= n
        n -= 1
    prod /= factorial(k)
    return prod

def main():
    parser = argparse.ArgumentParser(
        description='Calculate discrete convolutions of sum of three RVs,'
        ' each ranging from 1 to 6.')
    parser.add_argument("-n", action="store", dest="n", type=int, required=True)
    parser.add_argument("-pi", action="store", dest="pi", nargs=6,
                        metavar=('p1', 'p2', 'p3', 'p4', 'p5', 'p6'),
                        help="The probabilities of observing 1 through 6",
                        type=float, required=True)

    args = parser.parse_args()
    n = args.n
    pi = args.pi
    assert(sum(pi) == 1.0)

    h = h_Y(n, pi)
    err = 1.0 - sum(h)
    assert(err * err <= 10 ** -10)
    min10_values = min10(n, h)
    max10_values = max10(n, h)

    print "Min probabilities:"
    for y in min10_values:
        print (str(y+n) + " %.4g" % h[y])
    print "\nMax probabilities:"
    for y in max10_values:
        print (str(y+n) + " %.4g" % h[y])
    with open("h_%d_dp.csv" % n, "wb") as fil:
        for i in range(n, n + len(h)):
            s = str(i) + ",%.4e\n" % h[i-n]
            fil.write(s)
    mean, var = compute_stats(n, h)
    print "\nMean is %.9g" % mean
    print "Variance is %.4g" % var

    # Plot distributions
    # y = np.zeros(5*n+1)
    y = np.zeros([100])
    for i in range(len(y)):
        y[i] += i+n

    # Negative binomial approximation g
    phi = (mean/n) ** -1
    g = np.zeros(len(y))
    for i in range(len(y)):
        t1 = comb(y[i]-1, 49)
        t2 = (1-phi)**(y[i]-50)
        t3 = phi**50
        g[i] = t1*t2*t3

    # Normal distribution
    nvar = 0
    for i in range(len(pi)):
        nvar += pi[i] * ((i+1-mean/n)**2)
    nvar *= n
    print(nvar)
    n = mlab.normpdf(y, mean, math.sqrt(nvar))

    plt.title('Probability Distribution Approximations')
    # plt.plot(y, h, color='red', label='Our DP Distribution')
    plt.plot(y, h[:100], color='red', label='Our DP Distribution')
    plt.plot(y, g, color='blue', label='Negative Binomial')
    plt.plot(y, n, color='green', label='Normal')
    plt.legend()
    plt.show()


if __name__ == "__main__":
    main()
