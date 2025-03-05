'''
First Created on Sat May 9 18:34:00 2020 by
@author: thurcni@163.com, xzhang@tsinghua.edu.cn
Modified on Tue March 4 2025 by
@auther: liujk0725@outlook.com

Usage:
    from FEMModel import FEMModel
    model = FEMModel()
    model.load('path to model.json')
    model.plot()
    model.assemble()
    model.solve()
    model.print_stress()
    model.plot()
'''

import numpy as np
import json
from typing import List
import matplotlib.pyplot as plt

class FEMModel:
    '''
    Member variables defining the FEM model
    Title: (str) Title of the problem to be solved.
    nsd  : (int) Number of space dimensions.
    ndof : (int) Number of degrees-of-freedom per node.
    nnp  : (int) Number of nodal points.
    nel  : (int) Number of elements.
    nen  : (int) Number of element nodes.
    neq  : (int) Number of equations (D.O.F)
    nd   : (int) Number of nodes on the essential boundary.

    CArea: (numpy.array(nel)) Element values of cross-sectional area.
    E    : (numpy.array(nel)) Element values of Young's modulus.
    leng : (numpy.array(nel)) Element values of length
    stress:(numpy.array(nel)) Element values of stress

    x    : (numpy.array(nnp))x coordinate.
    y    : (numpy.array(nnp))y coordinates.
    IEN  : (numpy.array(nen,nel)) Element connectivity array.
    LM   : (numpy.array(nen,nel)) Location matrix.
    K    : (numpy.array(neq,neq)) Global stiffness matrix
    f    : (numpy.array(neq,1)) Global nodal force vector            
    d    : (numpy.array(neq,1)) Solution vector

    plot_truss: (bool) Plot truss ?
    plot_node : plot node number ?
    plot_tex  : Convert figures into PGFPlots figures in LaTex file ?
    '''
    Title = ""
    nsd   = 0
    ndof  = 0
    nnp   = 0
    nel   = 0
    nen   = 0
    neq   = 0
    nd    = 0

    CArea = np.array([])
    E     = np.array([])
    leng  = np.array([])
    stress= np.array([])

    x     = np.array([])
    y     = np.array([])
    IEN   = np.array([[]])
    LM    = np.array([[]])
    K     = np.array([[]])
    f     = np.array([[]])
    d     = np.array([[]])

    plot_truss = False
    plot_node  = False
    plot_tex   = False

    def __init__(self):
        pass
  
    def load(self, model_path: str)->None:
        """ 
        Load the model data from a json file.
        Args:
            model_path (str): Path to the json file.
        """

        with open(model_path) as f:
            model_data = json.load(f)

        self.Title = model_data['Title']
        self.nsd   = model_data['nsd']
        self.ndof  = model_data['ndof']
        assert self.nsd == self.ndof, "The number of space dimensions is not consistent with the number of degrees-of-freedom"
        self.nnp   = model_data['nnp']
        self.nel   = model_data['nel']
        self.nen   = model_data['nen']    
        self.neq   = self.ndof * self.nnp
        self.nd    = model_data['nd']

        # initialize K, d and f 
        self.f = np.zeros((self.neq, 1))            
        self.d = np.zeros((self.neq, 1))        
        self.K = np.zeros((self.neq, self.neq))

        # define the mesh
        self.x = np.array(model_data['x']).astype(np.float64)
        self.y = np.array(model_data['y']).astype(np.float64)
        if self.nsd == 3:
            self.z = np.array(model_data['z']).astype(np.float64)
        else:
            self.z = np.zeros((self.nnp,))
            # should not be used in 2D
        self.IEN = np.array(model_data['IEN'], dtype=np.int64)
        if self.IEN.shape[1] != self.nen:
            self.IEN = self.IEN.transpose()
        assert self.IEN.shape[1] == self.nen, "The number of nodes per element is not consistent with the input data"
        self.LM = np.zeros((self.nen * self.ndof, self.nel), dtype=np.int64)
        for e in range(self.nel):
            for j in range(self.nen):
                for m in range(self.ndof):
                    ind = j*self.ndof + m
                    self.LM[ind, e] = self.ndof*(self.IEN[e, j] - 1) + m
        # element and material data (given at the element)
        self.E     = np.array(model_data['E'])
        self.CArea = np.array(model_data['CArea'])
        self.leng  = np.sqrt(np.power(self.x[self.IEN[:, 1] - 1] - 
                                    self.x[self.IEN[:, 0] - 1], 2) +
                            np.power(self.y[self.IEN[:, 1] - 1] - 
                                    self.y[self.IEN[:, 0] - 1], 2))
        self.stress = np.zeros((self.nel,))

        # prescribed forces
        fdof = model_data['fdof']
        force = model_data['force']
        for ind, value in enumerate(fdof):
            self.f[value - 1][0] = force[ind]

        # output plots
        self.plot_truss = model_data['plot_truss']
        self.plot_node  = model_data['plot_node']
        # self.plot_tex   = model_data['plot_tex']
        self.plot_tex = False

    def plot(self, savefig: bool = False):
        '''
        Plot the truss structure.
        Args:
            savefig (bool): Save the figure to a file.
        '''
        fig = plt.figure()
        if self.plot_truss == "yes":
            if self.ndof == 1:
                x_dis = self.x + self.d.reshape(-1, self.ndof)[:, 0] * 1e4
                for i in range(self.nel):
                    XX = np.array([x_dis[self.IEN[i, 0]-1], 
                                x_dis[self.IEN[i, 1]-1]])
                    YY = np.array([0.0, 0.0])
                    plt.plot(XX, YY, "blue")

                    if self.plot_node == "yes":
                        plt.text(XX[0], YY[0], str(self.IEN[i, 0]))
                        plt.text(XX[1], YY[1], str(self.IEN[i, 1]))
            elif self.ndof == 2:
                x_dis = self.x + self.d.reshape(-1, self.ndof)[:, 0] * 1e4
                y_dis = self.y + self.d.reshape(-1, self.ndof)[:, 1] * 1e4
                for i in range(self.nel):
                    XX = np.array([x_dis[self.IEN[i, 0]-1], 
                                x_dis[self.IEN[i, 1]-1]])
                    YY = np.array([y_dis[self.IEN[i, 0]-1], 
                                y_dis[self.IEN[i, 1]-1]])
                    plt.plot(XX, YY, "blue", linewidth=1)

                    if self.plot_node == "yes":
                        plt.text(XX[0], YY[0], str(self.IEN[i, 0]))
                        plt.text(XX[1], YY[1], str(self.IEN[i, 1]))
            elif self.ndof == 3:
                from mpl_toolkits.mplot3d import Axes3D
                ax = fig.add_subplot(111, projection='3d')
                x_dis = self.x + self.d.reshape(-1, self.ndof)[:, 0] * 1e4
                y_dis = self.y + self.d.reshape(-1, self.ndof)[:, 1] * 1e4
                z_dis = self.z + self.d.reshape(-1, self.ndof)[:, 2] * 1e4 

                for i in range(self.nel):
                    XX = np.array([x_dis[self.IEN[i, 0]-1], 
                                x_dis[self.IEN[i, 1]-1]])
                    YY = np.array([y_dis[self.IEN[i, 0]-1], 
                                y_dis[self.IEN[i, 1]-1]])
                    ZZ = np.array([z_dis[self.IEN[i, 0]-1], 
                                z_dis[self.IEN[i, 1]-1]])
                    
                    ax.plot(XX, YY, ZZ, "blue", linewidth=1)

                    if self.plot_node == "yes":
                        ax.text(XX[0], YY[0], ZZ[0], str(self.IEN[i, 0]), color="red")
                        ax.text(XX[1], YY[1], ZZ[1], str(self.IEN[i, 1]), color="red")

                ax.set_xlabel("X Axis")
                ax.set_ylabel("Y Axis")
                ax.set_zlabel("Z Axis")
                ax.set_title("3D Truss Structure")

            else:
                raise ValueError("The dimension (ndof = {0}) given for the \
                                plottruss is invalid".format(self.ndof))
            
            plt.title("Truss Plot")
            plt.xlabel(r"$x$")
            plt.ylabel(r"$y$")
            if savefig:
                plt.savefig("truss_plot.png")
            plt.show()

            # Convert matplotlib figures into PGFPlots figures stored in a Tikz file, 
            # which can be added into your LaTex source code by "\input{fe_plot.tex}"
            if self.plot_tex == "yes":
                import tikzplotlib
                tikzplotlib.clean_figure()
                tikzplotlib.save("fe_plot.tex")
        
        print("\t2D Truss Params \n")
        print(self.Title + "\n")
        print("No. of Elements  {0}".format(self.nel))
        print("No. of Nodes     {0}".format(self.nnp))
        print("No. of Equations {0}".format(self.neq))

    def print_stress(self):
        '''
        Calculate and print stresses of every element
        '''
        for e in range(self.nel):
            de = self.d[self.LM[:, e]]  # nodal displacements for each element
            const = self.E[e]/self.leng[e]

            if self.ndof == 1:
                self.stress[e] = const*(np.array([-1, 1])@de)
            elif self.ndof == 2:
                IENe = self.IEN[e] - 1
                xe = self.x[IENe]
                ye = self.y[IENe]
                s = (ye[1] - ye[0])/self.leng[e]
                c = (xe[1] - xe[0])/self.leng[e]
                self.stress[e] = const*(np.array([-c, -s, c, s])@de)
            elif self.ndof == 3:
                IENe = self.IEN[e] - 1
                xe = self.x[IENe]
                ye = self.y[IENe]
                ze = self.z[IENe]
                x_cos = (xe[1] - xe[0])/self.leng[e]
                y_cos = (ye[1] - ye[0])/self.leng[e]
                z_cos = (ze[1] - ze[0])/self.leng[e]
                self.stress[e] = const*(np.array([-x_cos, -y_cos, -z_cos, x_cos, y_cos, z_cos])@de)
            else:
                raise ValueError("The dimension (ndof = {0}) given for the \
                                problem is invalid".format(self.ndof))
            print("{0}\t\t\t{1}".format(e+1, self.stress[e]))
        

    def assemble(self):
        '''
        Calculate and assemble element stiffness matrix
        '''
        for e in range(self.nel):
            const = self.CArea[e]*self.E[e]/self.leng[e]

            # calculate element stiffness matrix
            if self.ndof == 1:
                ke = const*np.array([[1, -1], [-1, 1]])
            elif self.ndof == 2:
                IENe = self.IEN[e] - 1
                xe = self.x[IENe]
                ye = self.y[IENe]
                s = (ye[1] - ye[0])/self.leng[e]
                c = (xe[1] - xe[0])/self.leng[e]

                s_s = s*s
                c_c = c*c
                c_s = c*s

                ke = const*np.array([[c_c, c_s, -c_c, -c_s],
                                    [c_s, s_s, -c_s, -s_s],
                                    [-c_c, -c_s, c_c, c_s],
                                    [-c_s, -s_s, c_s, s_s]])
            elif self.ndof == 3:
                IENe = self.IEN[e] - 1
                xe = self.x[IENe]
                ye = self.y[IENe]
                ze = self.z[IENe]
                delta_x = xe[1] - xe[0]
                delta_y = ye[1] - ye[0]
                delta_z = ze[1] - ze[0]
                T = np.array([[delta_x, delta_y, delta_z, 0, 0, 0],
                            [0, 0, 0, delta_x, delta_y, delta_z]])
                K = np.array([[1, -1], [-1, 1]])
                ke = const * T.T @ K @ T / (self.leng[e])**2
            else:
                raise ValueError("The dimension (ndof = {0}) given for the problem \
                                is invalid".format(self.ndof))
            self.K[np.ix_(self.LM[:,e], self.LM[:,e])] += ke
   
    def solve(self, method = "reduce"):
        '''
        Solve the system of equations
        '''
        assert method in ["reduce", "penalty"]
        if method == 'reduce':
            nd = self.nd
            neq = self.neq
            K_E = self.K[0:nd, 0:nd]
            K_F = self.K[nd:neq, nd:neq]
            K_EF = self.K[0:nd, nd:neq]
            f_F = self.f[nd:neq]
            d_E = self.d[0:nd]
            
            # solve for d_F
            d_F = np.linalg.solve(K_F, f_F - K_EF.T @ d_E) 

            # reconstruct the global displacement d
            self.d = np.append(d_E, d_F)
            
            # compute the reaction r
            f_E = K_E @ d_E + K_EF @ d_F
            
            # write to the workspace
            print('\nsolution d')
            print(self.d)
            print('\nreaction f =', f_E)
            
            return f_E
        if method == "penalty":
            penalty = np.diag(self.K).mean() * 1e7
            K = self.K.copy()
            for i in range(self.nd):
                K[i, i] = penalty
                self.f[i] = 0
            print(np.diag(K))
            self.d = np.linalg.solve(K, self.f)
            print('\nsolution d')
            print(self.d)
            f_E = (self.K @ self.d)[:self.nd]
            print('\nreaction f =', f_E)
            return f_E
