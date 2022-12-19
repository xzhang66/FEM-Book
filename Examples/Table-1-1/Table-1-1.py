#!/usr/bin/env python3
"""
Author : Xiong Zhang
Date   : 2022-06-28
Purpose: Calculate Pi by Inscribed Polygons and extrapolated by Wynn Epsilon method
"""

import numpy as np

def WynnEpsilon(sn, k):
    """ 
    Perform Wynn Epsilon Algorithm 
    """

    n = 2 * k + 1
    e = np.zeros((n+1, n+1))

    for i in range(1, n+1):
        e[i, 1] = sn[i-1]

    for i in range(3, n+2):
        for j in range(3, i+1):
            e[i-1, j-1] = e[i-2, j-3] + 1/(e[i-1, j-2] - e[i-2, j-2])

    ek = e[:, 1:n + 1:2]
    return ek

# --------------------------------------------------
def main():

    n  = np.logspace(0,8,9,base=2).astype(int)
    pn = n*np.sin(np.pi/n)

    pw = np.zeros(4)
    for i in range(1,5):
        en = WynnEpsilon(pn,i)
        pw[i-1] = en[-1, -1]

    print ("{:<5} {:<20} {:<20}".format('n','Pi-n','Pi-Wynn'))
    for k in range(n.size):
        if (k%2 == 0 and k>0):
            i = int(k/2)
            print("{:<5} {:.15f}    {:.15f}".format(n[k], pn[k], pw[i-1]))
        else:
            print ("{:<5} {:.15f}".format(n[k], pn[k]))

# --------------------------------------------------
if __name__ == '__main__':
    main()
    