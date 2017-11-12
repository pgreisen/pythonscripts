import numpy as np
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d.axes3d as axes3d


fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')


#data
#fname = ["acb","fgj","asd","tert","hjk","aer","cvb","ghj","qwf","hng","cde","iolk","xsa","ong","aasd","ghgh","rtrt","wewe","qwqw","bnbn","jhgd","aghty","dfg","derA","eryh","asdq","Dbgd","SsadI"]
#fx = [0.673574075,0.727952994,0.6745,0.793558,0.37664721,0.67939197,0.748490102,0.722048276,0.767189284,0.768296082,0.732383162,0.748429373,0.751570337,0.698977802,0.703398777,0.802607746,0.786242535,0.366661714,0.792490268,0.698636545,0.769904493,0.762656928,0.478595152,0.759151743,0.728607906,0.778099194,0.728575153,0.703794547]
#fy = [0.21,0.49,0.36,0.304,0.1009,0.329287287,0.752,0.309,0.277462605,0.333935925,0.326699919,0.72242944,0.358848707,0.298222369,0.03486,0.43058538,0.373973649,0.12288,0.9444,0.2494,0.384761779,0.382446444,0.35914706,0.298360515,0.391041147,0.363895412,0.312359532,0.343344197]
#fz = [18,8.620699842,9.2,17.40183,24.44101897,18,8.075948931,19,22.24192367,15,13.28436202,4.472128831,14.53195939,15.93922217,0,0,11.28194885,0,0,26.12423918,9.200498046,14.01392223,14.14545413,17.8320704,8.985897324,10.53443457,12.48561226,11.80438073]


##fname = ["FEN15","FEN14","FEN17","FEN16","FEN11","FEN10","FEN13","FEN12","FEN51","FEN50","FEN39","FEN38","FEN19","FEN18","FEN57","FEN56","FEN54","FEN4","FEN7","FEN30","FEN1","FEN3","FEN2","FEN37","FEN9","FEN8","FEN36","FEN31","FEN35","FEN61","FEN58","FEN34","FEN60","FEN42","FEN62","FEN32","FEN53","FEN43","FEN47","FEN40","FEN48","FEN41","FEN28","FEN29","FEN24","FEN25","FEN26","FEN27","FEN20","FEN21","FEN22","FEN23","FEN52","FEN44","FEN45","FEN5","FEN33","FEN55","FEN59","FEN6","FEN46","FEN49"]
# SC
##fx = [0.621329,0.72712,0.765821,0.679956,0.617988,0.599424,0.660564,0.709588,0.719612,0.652964,0.755742,0.610231,0.582499,0.665569,0.694257,0.666918,0.808658,0.591735,0.673421,0.729274,0.732103,0.65648,0.626306,0.670004,0.686148,0.71187,0.614255,0.72606,0.715682,0.663657,0.630007,0.676882,0.661256,0.568939,0.536255,0.732511,0.706636,0.624979,0.673088,0.656978,0.638881,0.557923,0.636926,0.673523,0.718163,0.647069,0.800019,0.738814,0.685791,0.716046,0.736019,0.771875,0.613889,0.734551,0.721603,0.596165,0.734966,0.759357,0.719416,0.643121,0.6152,0.66343]
##fy = [0.498151,0.633396,0.671416,0.568827,0.58175,0.548562,0.635933,0.614541,0.685178,0.79673,0.642902,0.701897,0.604532,0.626315,0.658236,0.66619,0.692121,0.464954,0.531101,0.622178,0.627424,0.586304,0.610228,0.681738,0.594352,0.558563,0.702772,0.632663,0.694695,0.717052,0.646735,0.668241,0.706944,0.624011,0.541268,0.70706,0.635004,0.634222,0.623198,0.120792,0.596891,0.585491,0.664694,0.633692,0.682404,0.654474,0.704637,0.684369,0.68373,0.661277,0.586887,0.597622,0.548812,0.684738,0.674011,0.538313,0.704685,0.60449,0.704396,0.603636,0.582191,0.81725]





# ddG
##fy = [-14.0523,-15.0291,-11.9993,-11.8569,-9.70304,-10.518,-16.3094,-16.7872,-17.0962,-16.0659,-14.4879,-18.1896,-18.308,-7.22715,-13.1411,-15.7073,-18.666,-13.8829,-13.8247,-15.3812,-16.9373,-14.1901,-14.2288,-16.5374,-14.0054,-9.4565,-15.9814,-13.7515,-15.2052,-14.3,-10.8235,-15.7438,-13.4416,-12.032,-12.0441,-15.4487,-10.861,-12.3312,-14.9253,-13.5446,-13.61,-3.96997,-13.4306,-13.6665,-15.169,-14.0697,-16.864,-15.4114,-9.90799,-15.7594,-14.7485,-15.0431,-12.6559,-11.7803,-14.032,-7.59315,-15.5399,-13.0045,-11.7771,-12.1474,-13.5243,-14.0821]

# sasa
##fy = [0.946848,0.92901,0.895593,0.958253,0.881752,0.974775,0.955063,0.960617,0.862774,0.864384,0.860705,0.938858,0.920843,0.867293,0.915493,0.949652,0.909708,0.960976,0.909412,0.897714,0.933941,0.936653,0.94886,0.901199,0.875072,0.988663,0.894203,0.877693,0.873419,0.905671,0.934949,0.922193,0.958169,0.900773,0.919818,0.866434,0.882105,0.920235,0.815818,0.86546,0.897752,0.926473,0.867429,0.882433,0.907616,0.95098,0.918484,0.845105,0.957944,0.953281,0.815288,0.880202,0.869992,0.826529,0.893792,0.939349,0.868783,0.852506,0.886153,0.905733,0.801708,0.932831]

# IFE
##fz = [-15.642,-14.9706,-16.7117,-19.1487,-13.8987,-17.1338,-16.3436,-16.6777,-17.1655,-16.0473,-15.4473,-18.5965,-17.565,-13.6175,-12.9061,-16.7276,-19.9318,-16.9116,-15.059,-16.6432,-19.2678,-19.1746,-16.4664,-16.9417,-17.4579,-17.5407,-16.3101,-14.9147,-15.2418,-14.9367,-12.8003,-15.9717,-16.9354,-13.0555,-11.7913,-15.9166,-10.8035,-13.7061,-15.8367,-13.8263,-13.5302,-11.4452,-15.212,-14.6361,-15.1816,-14.2992,-16.9616,-16.4491,-13.0772,-15.5934,-14.9511,-15.3306,-12.7058,-12.1146,-14.2251,-15.1667,-15.9262,-13.3428,-11.6218,-15.4775,-12.7161,-15.9142]


fname = ["FEN15","FEN14","FEN17","FEN16","FEN11","FEN10","FEN13","FEN12","FEN51","FEN50","FEN39","FEN38","FEN19","FEN18","FEN57","FEN56","FEN54","FEN4","FEN7","FEN30","FEN1","FEN3","FEN2","FEN37","FEN9","FEN8","FEN36","FEN31","FEN35","FEN61","FEN58","FEN34","FEN60","FEN42","FEN62","FEN32","FEN53","FEN43","FEN47","FEN48","FEN41","FEN28","FEN29","FEN24","FEN25","FEN26","FEN27","FEN20","FEN21","FEN22","FEN23","FEN52","FEN44","FEN45","FEN5","FEN33","FEN55","FEN59","FEN6","FEN46","FEN49"]
# SC
fx = [0.621329,0.72712,0.765821,0.679956,0.617988,0.599424,0.660564,0.709588,0.719612,0.652964,0.755742,0.610231,0.582499,0.665569,0.694257,0.666918,0.808658,0.591735,0.673421,0.729274,0.732103,0.65648,0.626306,0.670004,0.686148,0.71187,0.614255,0.72606,0.715682,0.663657,0.630007,0.676882,0.661256,0.568939,0.536255,0.732511,0.706636,0.624979,0.673088,0.638881,0.557923,0.636926,0.673523,0.718163,0.647069,0.800019,0.738814,0.685791,0.716046,0.736019,0.771875,0.613889,0.734551,0.721603,0.596165,0.734966,0.759357,0.719416,0.643121,0.6152,0.66343]
fy = [0.498151,0.633396,0.671416,0.568827,0.58175,0.548562,0.635933,0.614541,0.685178,0.79673,0.642902,0.701897,0.604532,0.626315,0.658236,0.66619,0.692121,0.464954,0.531101,0.622178,0.627424,0.586304,0.610228,0.681738,0.594352,0.558563,0.702772,0.632663,0.694695,0.717052,0.646735,0.668241,0.706944,0.624011,0.541268,0.70706,0.635004,0.634222,0.623198,0.596891,0.585491,0.664694,0.633692,0.682404,0.654474,0.704637,0.684369,0.68373,0.661277,0.586887,0.597622,0.548812,0.684738,0.674011,0.538313,0.704685,0.60449,0.704396,0.603636,0.582191,0.81725]
fz = [-15.642,-14.9706,-16.7117,-19.1487,-13.8987,-17.1338,-16.3436,-16.6777,-17.1655,-16.0473,-15.4473,-18.5965,-17.565,-13.6175,-12.9061,-16.7276,-19.9318,-16.9116,-15.059,-16.6432,-19.2678,-19.1746,-16.4664,-16.9417,-17.4579,-17.5407,-16.3101,-14.9147,-15.2418,-14.9367,-12.8003,-15.9717,-16.9354,-13.0555,-11.7913,-15.9166,-10.8035,-13.7061,-15.8367,-13.5302,-11.4452,-15.212,-14.6361,-15.1816,-14.2992,-16.9616,-16.4491,-13.0772,-15.5934,-14.9511,-15.3306,-12.7058,-12.1146,-14.2251,-15.1667,-15.9262,-13.3428,-11.6218,-15.4775,-12.7161,-15.9142]


#errorbar data
xSD = [0.041504064,0.02402152,0.059383144,0.038393713,0.054242278,0.018450667,0.083524242,0.042438697,0.036334793,0.023742101,0.041280224,0.003936522,0.025525758,0.031090602,0.027155833,0.038639074,0.061699064,0.11610088,0.075548578,0.059801071,0.069031082,0.071645685,0.050143938,0.049165738,0.020437116,0.046606225,0.039779165,0.019699934]
ySD = [0.015649804,0.12643117,0.068676131,0.016337,0.015050422,0.0651138,0,0.028590823,0.033705502,0.025962039,0,0,0.036646619,0.062000616,0,0,0.026584944,0.005923891,0,0.027485812,0,0.058142106,0.004978857,0.011233057,0.051596586,0.013837766,0,0.054340381]
zSD = [3.677693713,1.345712547,0.724095592,1.856309389,34.56482051,1.487978871,0,1.173906828,2.887602472,0.305603391,0,0,1.791653266,3.842020113,0,0,0.474818671,0,0,5.113750225,0,1.113374167,0.264111881,2.483847286,2.787214029,0.60047479,0,3.881040381]



#Class for 3d object
class thriidii:
    def __init__ (self, azimuut, elevation, x, y, z, d, n, gr, oy, oz, axesS, xl, yl, zl, prj, COL, randinp):
        
        
        self.AZ = azimuut 
        self.EL = elevation
    
        self.Dx = x
        self.Dy = y
        self.Dz = z
        
        #get limits
        
        self.maxDx = np.max(self.Dx)
        self.maxDy = np.max(self.Dy)
        self.maxDz = np.max(self.Dz)
        
        self.minDx = np.min(self.Dx)
        self.minDy = np.min(self.Dy)
        self.minDz= np.min(self.Dz)
        
        self.maxXYZ = np.max([np.max(self.Dx), np.max(self.Dy), np.max(self.Dz)])
        self.minXYZ = np.min([np.min(self.Dx), np.min(self.Dy), np.min(self.Dz)])        
        
        self.maxXY = np.max([np.max(self.Dx), np.max(self.Dy)])
        self.minXY = np.min([np.min(self.Dx), np.min(self.Dy)])        
        
        self.maxXZ = np.max([np.max(self.Dx), np.max(self.Dz)])
        self.minXZ = np.min([np.min(self.Dx), np.min(self.Dz)])        
        
        self.maxYZ = np.max([np.max(self.Dy), np.max(self.Dz)])
        self.minYZ = np.min([np.min(self.Dy), np.min(self.Dz)])
        
        print "MAX Dx", self.maxDx
        print "MAX Dy", self.maxDy 
        print "MAX Dz", self.maxDz
        
        print "MIN Dx", self.minDx 
        print "MIN Dy", self.minDy
        print "MIN Dz", self.minDz

        tmpsolution = ["FEN21", "FEN45", "FEN49"]
        
        for i in  np.arange(0, len(self.Dx)):
            

            if( n[i] in tmpsolution):
                ax.plot([self.Dx[i]], [self.Dy[i]], [self.Dz[i]], ls="None", marker="^", zorder=90, color="red", mec=COL)
            else:
                ax.plot([self.Dx[i]], [self.Dy[i]], [self.Dz[i]], ls="None", marker=".", zorder=90, color=COL, mec=COL)
            #plot 3d errorbars 
            ##ax.plot([self.Dx[i]-xSD[i], self.Dx[i]+xSD[i]], [self.Dy[i], self.Dy[i]], [self.Dz[i], self.Dz[i]], alpha=0.3, ls="-", marker="_", zorder=90, color=COL, mec=COL)
            ##ax.plot([self.Dx[i], self.Dx[i]], [self.Dy[i]-ySD[i], self.Dy[i]+ySD[i]], [self.Dz[i], self.Dz[i]], alpha=0.3, ls="-", marker="_", zorder=90, color=COL, mec=COL)
            ##ax.plot([self.Dx[i], self.Dx[i]], [self.Dy[i], self.Dy[i]], [self.Dz[i]-zSD[i], self.Dz[i]+zSD[i]], alpha=0.3, ls="-", marker="_", zorder=90, color=COL, mec=COL)
            
            
        #if gr = 1, plot names next to data
        ##PGif gr == 1:
        if 1 == 1:
            for i in np.arange(0, len(self.Dx)):
                # plot data points

                if(n[i] in tmpsolution):
                    #ax.text(self.Dx[i], self.Dy[i] + oy, self.Dz[i] + oz, "%s" % (n[i]), size=10, color=COL, zorder=90)
                    ax.text(self.Dx[i], self.Dy[i] + oy, self.Dz[i] + oz, "%s" % (n[i]), size=10, color="red", zorder=90)
               
        else:
            print gr
        

        #function to plot projections (2d plots)    
        def tuudii(arname, asim, elev, axesss): 
            for j in set(arname): #iterate trough unique elements in name list
                
                iii = -1
                temp_array = []
                try:
                    while 1:

                        tmpsolution = ["FEN21", "FEN45", "FEN49"]
                        # if(j in tmpsolution ):
                        iii = arname.index(j, iii+1) #find all maching indexes for names
                        print "match at", iii, j
                        if(j in tmpsolution):
                            temp_array.append(iii)
                            #ax.plot([self.Dx[i]], [self.Dy[i]], [self.minDz], ms="x", color="red", zorder=200)
                        
                            
                    
                except ValueError:
                    tmp_2dx = []
                    tmp_2dy = []
                    tmp_2dz = []
                    
                    tmp_2dx_cont = []                  
                    tmp_2dy_cont = []
                    tmp_2dz_cont = []
                    
                    if axesss == 3:
                        minX = self.minXYZ
                        minY = self.minXYZ
                        minZ =self.minXYZ
                        maxX =self.maxXYZ
                        maxY =self.maxXYZ
                        maxZ =self.maxXYZ
                    if axesss == 2:
                        minX = self.minXY
                        minY = self.minXY
                        minZ =self.minDz
                        maxX =self.maxXY
                        maxY =self.maxXY
                        maxZ =self.maxDz                
                    if axesss == 1:
                        minX = self.minXZ
                        minY = self.minDy
                        minZ =self.minXZ
                        maxX =self.maxXZ
                        maxY =self.maxDy
                        maxZ =self.maxXZ
                        
                    if axesss == 0:
                        minX = self.minDx
                        minY = self.minYZ
                        minZ =self.minYZ
                        maxX =self.maxDx
                        maxY =self.maxYZ
                        maxZ =self.maxYZ
                    if axesss == 4:
                        minX = self.minDx
                        minY = self.minDy
                        minZ =self.minDz
                        maxX =self.maxDx
                        maxY =self.maxDy
                        maxZ =self.maxDz
                    else:
                        pass
                        
                    for i in temp_array:
                    
                        tmp_2dx.append(self.Dx[i])
                        tmp_2dy.append(self.Dy[i])
                        tmp_2dz.append(self.Dz[i])
                        #depending of plotting angle choose where to plot projections
                        if asim < 90:
                            tmp_2dx_cont.append(minX)
                            x_cont = minX
                            tmp_2dy_cont.append(minY)
                            y_cont = minY
                        else:
                            tmp_2dx_cont.append(maxX)
                            x_cont = maxX
                            tmp_2dy_cont.append(minY)
                            y_cont = minY
                        
                        if elev > 0:
                            tmp_2dz_cont.append(minZ)
                            z_cont = minZ
                        else:
                            tmp_2dz_cont.append(maxZ)
                            z_cont = maxZ

                        
                    ax.plot(tmp_2dx, tmp_2dy, tmp_2dz_cont, ls="dotted", color="#C0C0C0")
                    
                    
                    ax.plot(tmp_2dx, tmp_2dy_cont, tmp_2dz, ls="dotted", color="#C0C0C0")
                    
                    
                    ax.plot(tmp_2dx_cont, tmp_2dy, tmp_2dz, ls="dotted", color="#C0C0C0")
                    
                    

                    for k in temp_array:
                        print "K:", k
                        ax.plot([x[k]], [y[k]], [z_cont], marker="o", color="#C0C0C0")
                        
                        ##ax.plot([x[k]-xSD[k], x[k]+xSD[k]], [y[k], y[k]], [z_cont, z_cont], alpha=0.3,marker="_", color="#C0C0C0")
                        ##ax.plot([x[k], x[k]], [y[k]-ySD[k], y[k]+ySD[k]], [z_cont, z_cont], alpha=0.3,marker="_", color="#C0C0C0")
                        
                        
                        
                        
                        ax.plot([x[k]], [y_cont], [z[k]], marker="o", color="#C0C0C0")
                        ##ax.plot([x[k]-xSD[k], x[k]+xSD[k]], [y_cont, y_cont], [z[k], z[k]], marker="_", color="#C0C0C0")
                        ##ax.plot([x[k], x[k]], [y_cont, y_cont], [z[k]-zSD[k], z[k]+zSD[k]], marker="_", color="#C0C0C0")
                        
                                         
                        
                        ax.plot([x_cont], [y[k]], [z[k]], marker="o", color="#C0C0C0")
                        
                        ##ax.plot([x_cont, x_cont], [y[k]-ySD[k], y[k]+ySD[k]], [z[k], z[k]], marker="_", color="#C0C0C0")
                        ##ax.plot([x_cont, x_cont], [y[k], y[k]], [z[k]-zSD[k], z[k]+zSD[k]], marker="_", color="#C0C0C0")



                        #ax.text(self.Dx[k], self.Dy[k]+0.01, z_cont+0.5,  "%s" % (n[k]), size=9, zorder=1, color="#C0C0C0")
                        #ax.text(self.Dx[k], y_cont+0.01, self.Dz[k]+0.5,  "%s" % (n[k]), size=9, zorder=1, color="#C0C0C0")
                        #ax.text(x_cont, self.Dy[k]+0.01, self.Dz[k]+0.5,  "%s" % (n[k]), size=9, zorder=1, color="#C0C0C0")
                        if(k == 48):
                            ax.text(self.Dx[k], self.Dy[k] - 0.01, z_cont + 0.5, "%s" % (n[k]), size=9, zorder=1,color="black")
                        else:
                            ax.text(self.Dx[k], self.Dy[k] + 0.02, z_cont + 0.5, "%s" % (n[k]), size=9, zorder=1,color="black")
                        ax.text(self.Dx[k], y_cont + 0.01, self.Dz[k] + 0.5, "%s" % (n[k]), size=9, zorder=1, color="black")
                        ax.text(x_cont, self.Dy[k] + 0.01, self.Dz[k] + 0.5, "%s" % (n[k]), size=9, zorder=1, color="black")

        ax.set_xlabel(xl)
        ax.set_ylabel(yl)
        ax.set_zlabel(zl)
        
        #uncomment next 2 lines to draw a viewing angle information on the plot
        
        #TITLE = "az: %s, el: %s" % (self.AZ, self.EL)
        #plt.title(TITLE)
            
        ax.azim = self.AZ
        ax.elev = self.EL
        


            
            
        if prj == 1:
            tuudii(n, self.AZ, self.EL, axesS)
            
            
        else:
            pass
        
        
        if axesS == 4:
            ax.set_xlim3d(self.minDx,self.maxDx)
            ax.set_ylim3d(self.minDy,self.maxDy)
            ax.set_zlim3d(self.minDz,self.maxDz)
            print "XYZ min:",self.minXYZ,"XYZ max:",self.maxXYZ
        if axesS == 3:
            ax.set_xlim3d(self.minXYZ,self.maxXYZ)
            ax.set_ylim3d(self.minXYZ,self.maxXYZ)
            ax.set_zlim3d(self.minXYZ,self.maxXYZ)
            print "XYZ min:",self.minXYZ,"XYZ max:",self.maxXYZ
        if axesS == 2:
            ax.set_xlim3d(self.minXY,self.maxXY)
            ax.set_ylim3d(self.minXY,self.maxXY)
            ax.set_zlim3d(self.minDz,self.maxDz)
            print "XY min:",self.minXYZ,"XY max:",self.maxXYZ, "Z min:", self.minDz,"Z max:", self.maxDz        
        if axesS == 1:
            ax.set_xlim3d(self.minXZ,self.maxXZ)
            ax.set_ylim3d(self.minDy,self.maxDy)
            ax.set_zlim3d(self.minXZ,self.maxXZ)
            print "XZ min:",self.minXZ,"XZmax:",self.maxXZ, "Y min:", self.minDy,"Y max:", self.maxDy 
        if axesS == 0:
            ax.set_xlim3d(self.minDx,self.maxDx)
            ax.set_ylim3d(self.minYZ,self.maxYZ)
            ax.set_zlim3d(self.minYZ,self.maxYZ)
            print "YZ min:",self.minYZ,"YZmax:",self.maxYZ, "X min:", self.minDx,"X max:", self.maxDx
        else:
            pass
        #plt.savefig("%s.png" % (randinp), format="png")
        plt.show()
        
# azimuut, elevation, x, y, z, ,name_list, gr, norm, Ylabeloffset, Zlabeloffset, xlabel, ylabel, zlabel, projections, color, random number)

#with projections
thriidii(45, 22, fx, fy, fz, fz, fname, 1,0.01, 0.5, 4, "SC", "pstat", "IFE", 1, "green",1)

#without projections
#thriidii(45, 22, fx, fy, fz, fz, fname, 1,0.01, 0.2, 0.2, "SC", "pstat", "IFE", 0, "green",1)
#init( azimuut, elevation, x, y, z, d, n, gr, oy, oz, axesS, xl, yl, zl, prj, COL, randinp):