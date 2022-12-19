#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Plot results for example 3-1 

Created on Wed Jun 24 17:47:07 2020

@author: xzhang
"""
import numpy as np
import tikzplotlib
import matplotlib.pyplot as plt

def exact(a, x):
    """ exact solution of Example 3-1"""
    b = np.sqrt(np.sqrt(a/4))
    division = a*(np.cos(b)**2*np.cosh(b)**2+np.sin(b)**2*np.sinh(b)**2)
    coef_A = np.cos(b)*np.cosh(b)/division
    coef_B = np.sin(b)*np.sinh(b)/division
    result = -1/a + coef_A*np.cos(b*x)*np.cosh(b*x) + coef_B*np.sin(b*x)*np.sinh(b*x)
    return result

def w(a1, x):
    """ 例3-1中的近似解w(x) """
    N1 = -(5-x**2)*(1-x**2)/24
    return a1*N1

def collocation(a, x):
    """ solution obtained by collocation method """
    a1 = 1/(1+5*a/24)
    return w(a1, x)

def subdomain(a, x):
    """ solution obtained by subdomain method """
    a1 = 1/(1+2*a/15)
    return w(a1, x)


def galerkin(a, x):
    """ solution obtained by Galerkin method """
    a1 = 1/(1+31*a/189)
    return w(a1, x)
 
def leastsquare(a, x):
    """ solution obtained by least square method """
    a1 = (1+2*a/15)/(1+4*a/15+62*a**2/2835)
    return w(a1, x)

def leastsquarecollocation(a, x):
    """ solution obtained by least square collocation method """
    a1 = 1/(1+44*a/243)
    return w(a1, x)


if __name__ == "__main__":
    plt.rcParams['font.sans-serif']=['SimSun'] #用来正常显示中文标签 
    plt.rcParams['axes.unicode_minus']=False #用来正常显示负号

    print("      中点(x=0)处的挠度和相对误差(%)")
    print("      a      精确解       配点法       子域法     伽辽金法    最小二乘法 最小二乘配点法")
    x = np.linspace(-1, 1, 25)
    for a in [1,10,100,1000]:
        yc = collocation(a,x)
        ys = subdomain(a,x)
        yg = galerkin(a,x)
        yl = leastsquare(a,x)
        ylc = leastsquarecollocation(a,x)
        ye = exact(a,x)

        fig = plt.figure()
        ax = fig.add_subplot()
        ax.plot(x, yc, color='red', label='配点法', marker=".", markevery=5)
        ax.plot(x, ys, color='blue', label='子域法', marker="+", markevery=6)
        ax.plot(x, yg, color='black', label='伽辽金法', marker="v", markevery=4)
        ax.plot(x, yl, color='pink', label='最小二乘法', marker="^", markevery=5)
        ax.plot(x, ylc, color='orange', label='最小二乘配点法', marker="*", markevery=5)
        ax.plot(x, ye, color='green', label='解析解', marker="o", markevery=4)

    #   ax.set_title(u'挠度曲线')
        ax.set_xlabel(r'$x$')
        ax.set_ylabel(r'$w(x)$')
        ax.legend()
    
        plt.show()
#        tikzplotlib.clean_figure()
#        tikzplotlib.save("example-3-1-%s-plot.tex"%a)

        with open("Elasti-Beam-%s.dat"%a, "a") as f:
            f.truncate(0)
            f.write("      x     collocation  subdomain    galerkin  leastsquare  LSCollocation    exact\n")
            np.savetxt(f, np.c_[x, yc, ys, yg, yl, ylc, ye], fmt='%1.4e')
        
        mc = collocation(a,0)
        ms = subdomain(a,0)
        mg = galerkin(a,0)
        ml = leastsquare(a,0)
        mlc = leastsquarecollocation(a,0)
        me = exact(a,0)

        ec = 100*(mc-me)/me
        es = 100*(ms-me)/me
        eg = 100*(mg-me)/me
        el = 100*(ml-me)/me
        elc = 100*(mlc-me)/me

        print("%6d%12.6f%12.6f%12.6f%12.6f%12.6f%12.6f" % (a, me, mc, ms, mg, ml, mlc))
        print("%18s   (%7.3f)   (%7.3f)   (%7.3f)   (%7.3f)   (%7.3f)" % (" ", ec, es, eg, el, elc))
