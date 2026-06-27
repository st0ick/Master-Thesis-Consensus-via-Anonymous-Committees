from gmpy2 import *
import numpy as np
from quorums.quorum_sizes import cod_committee_cutoff_probabilities, print_to_string, bin_pdf_left_tail_area


def optimal_quorum(N, L, S, b, startCommitteeSize=0, stepSize=1000, measure=lambda x,y: max(x,y), print_output=False):
    target_prob = mpfr('2')**(-b)
    expectedCommitteeSize = startCommitteeSize
    p = mpfr('1')
    #Find smallest quorum Q that satisfies the probability conditions.
    while(1):
        p = mpfr(expectedCommitteeSize / N) 
        
        Q, r_tail, l_tail = cod_committee_cutoff_probabilities(L, S, p)
        #temp
        logL = gmpy2.log2(r_tail)
        logR = gmpy2.log2(l_tail)
        #
        if measure(r_tail, l_tail) < target_prob:
            if stepSize == 1:
                break 
            else:
                expectedCommitteeSize -= stepSize
                stepSize = stepSize // 2
        
        expectedCommitteeSize += stepSize 
    counter = 1
    expectedCommitteeSize_N = expectedCommitteeSize - 2

    # Check that optimal quorums for larger committee sizes also satisfy the probability constraints. Because of the discrete nature,
    # it is possible that some larger committee size leads to slightly larger probabilities.

    fileName = "Quorums.txt"

    f = open(fileName, "w")

    while(counter < 20):
        p_N = mpfr('1') * expectedCommitteeSize_N / N
        QN, r_tail_N, l_tail_N = cod_committee_cutoff_probabilities(L, S, p_N)

        # Print 
        if print_output:
            logL = gmpy2.log2(r_tail_N)
            logR = gmpy2.log2(l_tail_N)

            f.write(print_to_string("p = ", p_N, " (expected committee size = ", expectedCommitteeSize_N, ")", sep=''))
            f.write("  Q = %d\n" % (QN))
            f.write(print_to_string("  log(L) = ", logL, sep=''))
            f.write(print_to_string("  log(R) = ", logR, sep=''))
            f.write("\n")


        if r_tail_N < 0 or measure(r_tail_N, l_tail_N) >= target_prob:
            counter = 0
            Q, p, r_tail, l_tail = QN, p_N, r_tail_N, l_tail_N
            expectedCommitteeSize = expectedCommitteeSize_N + 1
        else:
            counter += 1

        expectedCommitteeSize_N += 1
        
    print(f'Optimal Q = {Q}, with expected committeesize = {expectedCommitteeSize}')
    return Q, p

def get_tail_areas(N, fM, fI, k, p, cod=False ):
    H = N * (1-fM-fI)
    M = N * fM
    T = M + H // 2

    L = H if cod else N*(1-fI) 
    
    r_tail = bin_pdf_left_tail_area(L, p, k-1)
    l_tail = 1 - bin_pdf_left_tail_area(T, p, k-1) 

    return r_tail, l_tail


# Find 
def func_fM(N, fM, fI, b, n_steps=35, cod=False, measure=lambda x,y: max(x,y)):
    H = N * (1-fM-fI)
    M = N * fM
    T = M + H // 2

    L = H if cod else N*(1-fI) 
    k, p = optimal_quorum(N, L, T, b, startCommitteeSize=1000,  measure=measure)

    r = np.linspace(0, 0.15, n_steps)


    x = [i for i in r]
    y = [get_tail_areas(N, fM + delta, fI, k, p, cod=cod) for delta in r]

    return x, y, k, p

def find_cutoff_fM(config, b_start, b_cutoff, fM_start, cod=False, measure=lambda x, y: max(x,y), input=-1):
    N, fM, fI = config 
    N = mpz(N)
    fM = mpfr(fM)
    fI = mpfr(fI)
    H = N * (1-fM-fI)
    M = N * fM
    T = M + H // 2
    L = H if cod else N*(1-fI) 

    Q, p = optimal_quorum(N, L, T, b_start) if input == -1 else input

    fM_start = fM + fM_start
    M_t = gmpy2.floor(fM_start * N)
    H_t = gmpy2.floor(N*(1-fM_start-fI))
    T_t = M_t + H_t // 2
    stepSize=1000000
    while(H_t > 0):
        fM_t = M_t / N 
        left_tail, right_tail = get_tail_areas(N, fM_t, fI, Q, p, cod)
        logL = gmpy2.log2(left_tail)
        logF = gmpy2.log2(right_tail)
        if measure(left_tail, right_tail) >= 2**(-b_cutoff):
            M_t -= stepSize 
            H_t += stepSize
            stepSize = stepSize // 2
            if stepSize == 1:
                fM_cutoff_l = (M_t - 1) / N
                fM_cutoff_r = M_t / N
                l_tail, r_tail = get_tail_areas(N, fM_cutoff_l, fI, Q, p, cod)
                return fM_cutoff_l - fM, fM_cutoff_r - fM, -gmpy2.log2(measure(l_tail, r_tail))
            stepSize = stepSize // 2
        else:
            M_t += stepSize 
            H_t -= stepSize 
        



def func_fI(N, fM, fI, b, n_steps=31, cod=False, measure=lambda x,y: max(x,y)):
    I = gmpy2.rint_floor(N * fI)
    
    M = gmpy2.rint_floor(N * fM)
    H = N - I - M
    T = M + H // 2
    L = H if cod else (N - I)
    
    k, p = optimal_quorum(N, L, T, b, measure=measure)

    r = np.linspace(-0.15, 0.15, n_steps)
    x = [i for i in r]
    y = []
    for delta in r:
        y.append(get_tail_areas(N, fM, fI + delta, k, p, cod))
    return x, y, k, p


def find_cutoff_fI_left(config, b_start, b_cutoff, dfI_start, cod=False, measure=lambda x, y: max(x,y), input=-1):
    N, fM, fI = config 
    N = mpz(N)
    fM = mpfr(fM)
    fI = mpfr(fI)
    H = N * (1-fM-fI)
    M = N * fM
    T = M + H // 2
    L = H if cod else N*(1-fI) 

    Q, p = optimal_quorum(N, L, T, b_start) if input == -1 else input

    fI_start = fI + dfI_start
    
    I_t = N*fI_start
    H_t = N - I_t - M
    stepSize = 1000000
    while H_t < N:
        fI_t = I_t / N 
        left_tail, right_tail = get_tail_areas(N, fM, fI_t, Q, p, cod)
        if measure(left_tail, right_tail) >= 2**(-b_cutoff):
            H_t -= stepSize 
            I_t += stepSize 
        
        else:
            H_t += stepSize
            I_t -= stepSize 
            
            if stepSize == 1:
                fI_cutoff_l = I_t / N
                fI_cutoff_r = (I_t + 1) / N
                l_tail, r_tail = get_tail_areas(N, fM, fI_cutoff_r, Q, p, cod)
                return fI_cutoff_l - fI, fI_cutoff_r - fI, -gmpy2.log2(measure(l_tail, r_tail))
            stepSize = stepSize // 2



   

def find_cutoff_fI_right(config, b_start, b_cutoff, dfI_start, cod=False, measure=lambda x, y: max(x,y), input=-1):
    N, fM, fI = config 
    N = mpz(N)
    fM = mpfr(fM)
    fI = mpfr(fI)
    H = N * (1-fM-fI)
    M = N * fM
    T = M + H // 2
    L = H if cod else N*(1-fI) 

    Q, p = optimal_quorum(N, L, T, b_start) if input == -1 else input

    fI_start = fI + dfI_start
    
    I_t = N*fI_start
    H_t = N - I_t - M
    stepSize = 1000000
    while I_t < N:
        fI_t = I_t / N 
        left_tail, right_tail = get_tail_areas(N, fM, fI_t, Q, p, cod)
        if measure(left_tail, right_tail) < 2**(-b_cutoff):
            H_t -= stepSize 
            I_t += stepSize 
        
        else:
            H_t += stepSize
            I_t -= stepSize 
            
            if stepSize == 1:
                fI_cutoff_l = I_t  / N
                fI_cutoff_r = (I_t + 1) / N
                l_tail, r_tail = get_tail_areas(N, fM, fI_cutoff_l, Q, p, cod)
                return fI_cutoff_l - fI, fI_cutoff_r - fI, -gmpy2.log2(measure(l_tail, r_tail))
            stepSize = stepSize // 2




