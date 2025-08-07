import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as ani
from celluloid import Camera
from os import getcwd
import time

def article_params():
    plt.rcParams.update({"text.usetex":True,"figure.figsize":[3.5,3],"font.family": "serif","font.sans-serif": ["Helvetica"],'font.size':12,'lines.linewidth':2,'legend.fontsize':11,'xtick.labelsize':10,'ytick.labelsize':10})
def viewing_params():
    plt.rcParams.update({"text.usetex":False,"figure.figsize":[12,9],"font.family": "serif","font.sans-serif": ["Helvetica"],'font.size':22,'lines.linewidth':3,'legend.fontsize':20,'xtick.labelsize':20,'ytick.labelsize':20})
viewing_params()
tabcolors = ('tab:blue','tab:orange','tab:green','tab:red','tab:purple','tab:brown','tab:pink','tab:gray','tab:olive','tab:cyan')
mpllinestyles = ('solid','dashed','dotted','dashdot',(0, (3, 5, 1, 5, 1, 5)))
mpllinestyles = ('solid','dashed','dotted','dashdot')
mplmarkers = ('o','v','^','<','>','s','*','D','P','X')

variable_dict = {'c':'concentration', 'c1':'left concentration', 'c2':'right concentration', 'u':'velocity', 'q':'q, conserved velocity', 'phi':'phi, conserved concentration', 'h':'height','h_latex':'$h(x,t)$','u_latex':'$u(x,t)$','c1_latex':'$c_1(x,t)$','c2_latex':'$c_2(x,t)$'}
char_to_cons = {'u':'q','c1':'phi1','c2':'phi2'}
char_vars = ['u','c1','c2']

class LoadSim:

    def __init__(self, hR0, cR0, U_s, rootFile='data/',VARS=['h','u','c1','c2'], subFile = ''):
        self.hR0 = hR0
        self.cR0 = cR0
        self.U_s = U_s
        self.rootFile = rootFile
        self.subFile = subFile
        
        self.fileName = "hTwo%0.2f_cTwo%0.2f_Us%0.3f"%(self.hR0,self.cR0,self.U_s)
        
        self.unpack(VARS)
        self.dt = self.T[1]-self.T[0]
        with open(self.rootFile + self.subFile + self.fileName + '/info.log') as fh:
            for line in fh:
                if 'Run time was ' in line: 
                    self.runTime = float(line[line.find('was')+4:line.find('seconds')-1])
                if 'collision time:' in line: 
                    self.coll_time = float(line[line.find('=')+1:])
                if 'collision index:' in line: 
                    self.coll_idx = int(float(line[line.find('=')+1:]))
                if 'collision position:' in line: 
                    self.coll_loc = float(line[line.find('=')+1:])

    def sim_info(self):
        # print the log file to screen. 
        print('')
        with open(self.rootFile + self.subFile + self.fileName + '/info.log') as fh:
            for line in fh: print(line)

    def unpack(self, whichVars):
        def single_unpack(whichFile):
            """
            Argument is the file name of the data
            When loading data, Each row is a ''time slice'' which is a vector of the function values at each point in space
        
            The first row is the initial time and then spatial vector of the center of the cells
            The first column is the time values, with the initial time listed twice, see above comment
            Removing the first column and first row you are left with a matrix of the function values (rows are time slices)
            """
            if whichFile in char_vars:
                A = np.loadtxt(self.rootFile + self.subFile + self.fileName + '/' + char_to_cons[whichFile])
            else:
                A = np.loadtxt(self.rootFile + self.subFile + self.fileName + '/' + whichFile)
            self.x = A[0,1:]
            self.T = A[1:,0]
            DependentVariable = A[1:,1:]
            if whichFile in char_vars:
                try: 
                    DependentVariable = DependentVariable/self.h
                except AttributeError:
                    self.h = single_unpack('h')
                    DependentVariable = DependentVariable/self.h
            return DependentVariable

        for var in whichVars:
            setattr(self,var,single_unpack(var))

    def makeMP4(self, varList = ['h','u','c1','c2'],tMax=1000.,show_legend = True, xlim = None, ylim = None,framerate = 30.):

        fig = plt.figure(figsize=(12,8))
        camera = Camera(fig)
    
        U = []
        for var in varList:
            U.append(getattr(self,var))

        if ylim:
            plt.ylim(ylim)
            ymin1,ymax1 = ylim
            ymin2,ymax2 = ylim
        else:
            ymin1,ymax1 = min(np.min(U[0]),np.min(U[1])),max(np.max(U[0]),np.max(U[1]))
            ymin2,ymax2 = min(np.min(U[2]),np.min(U[3])),max(np.max(U[2]),np.max(U[3]))
        maxVEL = 0.

        for ii in range(len(self.T)):
            if self.T[ii]>tMax: continue
            subplotcounter = 1
            for jj in range(len(varList)):
                subplotcounter += 1
                u = U[jj]
                plt.subplot(2,1,int(subplotcounter/2))
                plot_ = plt.plot(self.x,u[ii,:],color=tabcolors[jj],linestyle = mpllinestyles[0])
                plt.grid()
    
            plt.subplot(211)
            plt.xlim(xlim)
            plt.ylim([ymin1,ymax1])
            VEL = U[1]
            maxVEL = max(maxVEL,np.max(VEL[ii,:]))
            timeStr = 't = %0.1f'%(self.T[ii])
            plt.text((xlim[0]+xlim[1])/2 if xlim else 0., ymax1 + 0.1*(ymax1-ymin1),timeStr,verticalalignment = 'center',horizontalalignment = 'center')

            plt.subplot(212)
            plt.xlim(xlim)
            if np.abs(self.T[ii]-int(self.T[ii]))<self.dt/2: print('t = %i'%np.round(self.T[ii]))
            camera.snap()
    
        if show_legend:
            Legend = []
            for var in varList:
                Legend.append(variable_dict[var])
            plt.subplot(211)
            plt.legend(Legend[:2],loc= 'upper left')
            plt.subplot(212)
            plt.legend(Legend[2:],loc= 'upper left')
    
        animat = camera.animate()

        Writ = ani.FFMpegWriter(fps=framerate, metadata=dict(artist='nathan'))
        animat.save(getcwd() + '/' + self.fileName.replace('.','_') + '.mp4', writer = Writ)
        plt.close()

    def plot_time(self,var,desired_time):
        index = np.argmin(np.abs(self.T-desired_time))
        plt.plot(self.x,getattr(self,var)[index,:],label = '$t=%0.1f$'%self.T[index])

    def plot_times(self,var,times,xlim=None,show_legend = True):
        for t in times:
            self.plot_time(var,t)
        plt.xlim(xlim)
        plt.xlabel('x')
        plt.ylabel(variable_dict[var])
        if show_legend: plt.legend()
        plt.savefig(getcwd() + '/' + self.fileName.replace('.','_') + '.png')
        plt.show()

    def plot_examples(self,var,times=[0,4,8,12]):
        article_params()
        plt.figure(figsize = [8,4])
        for i,v in enumerate(var):
            plt.subplot(len(var),1,i+1)
            self.plot_times(v,times,xlim=[-10,10],wl=False)
            plt.ylabel(variable_dict[v + '_latex'])
            if i < len(var)-1: plt.gca().set_xticks([])
        
        plt.subplot(len(var),1,1)
        plt.legend(['$t=%0.1f$'%t for t in times],ncol=len(times),loc='upper center', bbox_to_anchor=(0.5, 1.6))
        plt.subplot(len(var),1,len(var))
        plt.xlabel('$x$')
        plt.savefig(self.rootFile + 'solutions/plots/solution_example_' + self.fileName + '.png', bbox_inches='tight',dpi=400)

    def deposition_details(self):
        dx = self.x[1]-self.x[0]

        self.xr = self.x[self.coll_idx+1:int(self.N*0.8)]
        self.l2r = self.d1[-1,self.coll_idx+1:int(self.N*0.8)]

        self.xl = self.x[int(self.N*0.2):self.coll_idx+1]
        self.r2l = self.d2[-1,int(self.N*0.2):self.coll_idx+1]

        two_sided_intrusion = np.hstack([self.r2l,self.l2r])
        two_sided_x = np.hstack([self.xl,self.xr])
        
        l2r_mass = np.sum(self.l2r)*dx
        r2l_mass = np.sum(self.r2l)*dx
 
        self.intrusion_mass = l2r_mass - r2l_mass
        self.COM_x = np.sum(two_sided_x*two_sided_intrusion)*dx/(l2r_mass + r2l_mass)

    def plot_deposit(self, LS = 'solid'):
        plt.plot(self.x,self.d1[-1,:],label='$d_1(x), U_s = %0.3f$'%self.U_s,color = 'tab:blue', linestyle = LS)
        plt.plot(self.x,self.d2[-1,:],label='$d_2(x), U_s = %0.3f$'%self.U_s,color = 'tab:orange', linestyle = LS)
        plt.xlim([-20,20])
        plt.title('$h_{2,0}$ = %0.2f, $c_{2,0}$ = %0.2f'%(self.hR0,self.cR0))
        plt.xlabel('$x$')

    def plot_intrusion(self,wantLegend=True, want_y_label=True):
        try:
            self.intrustion_mass
        except AttributeError:
            self.deposition_details()
        plt.plot(self.xl,self.r2l,label = 'right-to-left: $d_2(x), x>x_c$')
        plt.plot(self.xr,self.l2r,label = 'left-to-right: $d_1(x), x<x_c$')
        plt.plot([self.coll_loc]*2,[0,max(np.max(self.l2r),np.max(self.r2l))],color = 'k',linestyle = 'dashed',label = 'collision location $(x_c)$: %0.2f'%self.coll_loc,linewidth = 2)
        plt.plot([self.COM_x]*2,[0,max(np.max(self.l2r),np.max(self.r2l))],color = 'red',linestyle = 'dashed',label = 'COM $x$-coordinate: %0.2f'%self.COM_x)
        try: 
            xMin = self.xl[np.argwhere(self.r2l>0.01*np.max(self.r2l))[0][0]]
            xMax = self.xr[np.argwhere(self.l2r>0.01*np.max(self.l2r))[-1][0]]
            plt.xlim([xMin,xMax])
        except IndexError:
            print('Adaptive xlim did not work, setting bounds at [-5,5]')
            plt.xlim([-5,5])
        plt.title('$h_{2,0}$ = %0.2f, $c_{2,0}$ = %0.2f \n $\int_{x_c}^\infty d_1(x) dx - \int_{-\infty}^{x_c} d_2(x) dx $= %0.2f'%(self.hR0,self.cR0,self.intrusion_mass))
        if wantLegend: plt.legend()
        plt.xlabel('$x$')
        if want_y_label: plt.ylabel('deposits \n left: $d_1(x)$, right: $d_2(x)$')

    def sus_conc(self,LS = None, LC = None):
        c = (self.c1 + self.c2)*self.h
        dx = self.x[1]-self.x[0]
        initial = np.sum(c[0,1:] + c[0,:-1])*(dx/2)
        suspended_pct = np.sum(c[:,1:] + c[:,:-1],1)*(dx/2)/initial
        plt.plot(self.T,suspended_pct, color = LC, linestyle = LS)
        return self.T[np.argwhere(suspended_pct > 0.05)[-1][0]]

    def bore_data(self,subSampleBy=1,all_var=True,t_start=0,plot=False,save = True):
        try: 
            C=np.loadtxt(self.rootFile + 'solutions/postData/post_collision_bore_data_'+ self.fileName + '.csv',delimiter=',')
            t_post = C[:,0]
            bore = C[:,1]
            h_plus,   h_minus   = C[:,2], C[:,3]
            u_plus,   u_minus   = C[:,4], C[:,5]
            front = C[:,6]
        except OSError: 
            print('cannot open, so I will post process the data MYSELF!')
            t_post = []
            dt = self.T[1]-self.T[0]
            post_collision = False
            h_max = 0
            u_plus, u_minus = [],[]
            h_plus, h_minus = [],[]
            phi_plus, phi_minus = [],[]
            bore, front = [],[]
            bore_index=0
            THRESH = []
            t_count = 0
            x_ = self.x[int(self.N/2):]
            for t,u_fake,h_fake in zip(self.T,self.u[:,int(self.N/2):],self.h[:,int(self.N/2):]):
                u_ = deepcopy(u_fake)
                h_ = deepcopy(h_fake)
                window_size = 50
                if (not post_collision) and h_[0]>2*self.h_min and np.max(h_)<h_max and x_[np.argmax(h_)]<2 and t>t_start:
                    post_collision = True
                front_index = np.argwhere(h_>2*self.h_min)[-1][0]
                if post_collision:
                    threshold = (h_[int((front_index+bore_index)/2)] + h_[int(bore_index/2)])/2
                    THRESH.append(threshold)
                    u_max = np.max(u_)
                    if bore_index:
                         h_[:max(bore_index-2*window_size,0)]=np.max(h_)
                         h_[min(bore_index+2*window_size,len(h_)):]=self.h_min
                    bore_index = np.argwhere(h_>threshold)[-1][0]
                    u_bore = u_[max(bore_index-window_size,0):bore_index+window_size]
                    h_bore = h_[max(bore_index-window_size,0):bore_index+window_size]
                    x_bore = x_[max(bore_index-window_size,0):bore_index+window_size]
                    u_mbi, u_pbi = np.argmax(u_bore), np.argmin(u_bore)

                    bore_loc = (x_[bore_index]-x_[bore_index+1])*(threshold-h_[bore_index+1])/(h_[bore_index]-h_[bore_index+1])+x_[bore_index+1]
                    if bore_loc < 0.2: continue
                    if plot and t<2.5:
                        plt.subplot(211)
                        p1 = plt.plot(x_bore,h_bore,label='t=%0.2f'%t)[0]
                        h_boxed = np.array([h_bore[u_mbi] if x_b<bore_loc else h_bore[u_pbi] for x_b in x_bore])
                        plt.plot(x_bore,h_boxed,linestyle = 'dashed',color = p1.get_color())

                        plt.subplot(212)
                        plt.plot(x_bore,u_bore)
                     
                    t_post.append(t-self.coll_time)
                    u_minus.append(u_bore[u_mbi])
                    u_plus.append(u_bore[u_pbi])
                    h_minus.append(h_bore[u_mbi])
                    h_plus.append(h_bore[u_pbi])
                    bore.append(bore_loc)
                    front_loc = x_[front_index]
                    front.append(front_loc)
                if plot: 
                    plt.subplot(212)
                    plt.xlabel('x')
                    plt.ylabel('velocity')

                    plt.subplot(211)
                    plt.xlabel('x')
                    plt.ylabel('height')

                    plt.legend()
                h_max = np.max(h_)
    
            array_to_save=np.array((np.array(t_post),np.array(bore),np.array(h_plus),np.array(h_minus),np.array(u_plus),np.array(u_minus),np.array(front))).T
            if save: np.savetxt(self.rootFile + 'solutions/postData/post_collision_bore_data_' + self.fileName + '.csv',array_to_save,delimiter=',')
        def subsample(x,subSampleBy = subSampleBy):
            x = np.array(x)
            return x[range(0,x.shape[0],subSampleBy)]
        self.t_post = subsample(t_post)
        self.bore = subsample(bore)
        self.hP_data, self.hM_data = subsample(h_plus),subsample(h_minus)
        self.uP_data, self.uM_data = subsample(u_plus),subsample(u_minus)
        self.front_data = subsample(front)


class DepositionAnalysis:
    def __init__(self, U_s, rootFile, H2=np.linspace(0.7,1.42,73), C2=np.linspace(0.7,1.42,73), N = 5000, NuRe = 1000, finalTime = 40., h_min = 0.0001, CFL = 0.1, sharp = 50, apart = 5., FrSquared = 1., hL0 = 1.0, cL0 = 1.0, NuPe = None, subFile = 'solutions/postData/'):
        self.H2 = H2
        self.C2 = C2
        self.U_s = U_s
        self.rootFile = rootFile
        self.subFile = subFile
        self.N = N
        self.NuRe = NuRe
        self.finalTime = finalTime
        self.h_min = h_min
        self.CFL = CFL
        self.sharp = sharp
        self.apart = apart
        self.FrSquared = FrSquared
        self.hL0, self.cL0 = hL0, cL0
        
        self.fileName = "%iby%i_%iapart_N%i_CFL%0.3f_T%0.1f_NuRe%i_FrFr%0.3f_Us%0.3f_hmin%0.5f_sharp%i"%(self.H2.shape[0],self.C2.shape[0],self.apart,self.N,self.CFL,self.finalTime,self.NuRe,self.FrSquared,self.U_s,self.h_min,self.sharp)

        try: 
            self.intrusion_mass = np.loadtxt(self.rootFile + self.subFile + 'intrMass_' + self.fileName + '.csv', delimiter = ',')
            self.COM_x = np.loadtxt(self.rootFile + self.subFile + 'COMx_' + self.fileName + '.csv', delimiter = ',')
        except OSError:
            print('Cannot find processed sedimentation data. Post processing now.') 
            start = time.time()
            self.intrusion_mass = np.zeros([C2.shape[0],H2.shape[0]])
            self.COM_x = np.zeros([C2.shape[0],H2.shape[0]])
            self.deposited_mass = np.zeros([C2.shape[0],H2.shape[0]])
            self.coll_time = np.zeros([C2.shape[0],H2.shape[0]])
            self.coll_loc = np.zeros([C2.shape[0],H2.shape[0]])
            for i,h2 in enumerate(H2):
                for j,c2 in enumerate(C2):
                    if h2*c2<=1:
                        temp_sim = LoadSim(h2,c2,self.U_s,self.rootFile,['d1','d2'],N=self.N, NuRe=self.NuRe, finalTime=self.finalTime, h_min=self.h_min, CFL=self.CFL, sharp=self.sharp, apart=self.apart, FrSquared=self.FrSquared)
                        if temp_sim.d1.shape[0]>1 and temp_sim.d2.shape[0]>1:
                            temp_sim.deposition_details()
                            self.intrusion_mass[j,i] = temp_sim.intrusion_mass
                            self.COM_x[j,i] = temp_sim.COM_x
                            self.deposited_mass[j,i] = (np.sum(temp_sim.d1)+np.sum(temp_sim.d2))*(temp_sim.x[1]-temp_sim.x[0])/(temp_sim.hL0*temp_sim.cL0 + h2*c2)
                            self.coll_time[j,i] = temp_sim.coll_time
                            self.coll_loc[j,i] = temp_sim.coll_loc
                        else:
                            print('h2 = %0.2f, c2 = %0.2f, U_s = %0.3f simulation did not finish. Final deposit data DNE'%(h2,c2,self.U_s))
                            self.intrusion_mass[j,i] = np.nan
                            self.coll_time[j,i] = np.nan
                            self.coll_loc[j,i] = np.nan
                            self.COM_x[j,i] = np.nan
                    else:
                        self.intrusion_mass[j,i] = np.nan
                        self.COM_x[j,i] = np.nan
                        self.coll_time[j,i] = np.nan
                        self.coll_loc[j,i] = np.nan
            np.savetxt(self.rootFile + self.subFile + 'intrMass_' + self.fileName + '.csv', self.intrusion_mass, delimiter = ',')
            np.savetxt(self.rootFile + self.subFile + 'COMx_' + self.fileName + '.csv', self.COM_x, delimiter = ',')
            print('%0.2f seconds'%(time.time()-start))

    def myPcolor(self,attr,plotTitle,streamlines=True):
        pplot = plt.pcolormesh(self.H2,self.C2,getattr(self,attr),shading = 'gouraud')
        plt.contour(self.H2,self.C2,getattr(self,attr),colors='k',levels=[0])
        plt.colorbar(pplot)
        plt.gca().set_aspect('equal')
        plt.title(plotTitle)
        plt.xlabel('$h_{2,0}$')
        plt.ylabel('$c_{2,0}$')

    def plot_dimensional_analysis(self):
        article_params()
        plt.figure(figsize=[7,2.75])
        plt.subplot(121)
        self.myPcolor('intrusion_mass','Intrusion Mass')
        plt.subplot(122)
        self.myPcolor('COM_x','COM $x$-coordinate')
        plt.tight_layout()

        plt.savefig(self.rootFile + 'solutions/plots/deposition_analysis_' + self.fileName + '.png', bbox_inches='tight',dpi=400)
        plt.close()
        plt.rcParams.update({"text.usetex":False})

    def linear_appr(self,X,Y,Z):
        Xf,Yf,Zf = X.flatten(), Y.flatten(), Z.flatten()
        no_nans = ~np.isnan(Zf)
        Xf,Yf,Zf = Xf[no_nans],Yf[no_nans],Zf[no_nans]
        A = np.array([Xf,Yf,np.ones(Zf.shape[0])]).T
        a, b, c =  np.linalg.solve(np.matmul(A.T,A),np.matmul(A.T,np.array([Zf]).T)).flatten()
        print('a = %0.4f, b = %0.4f, c = %0.4f'%(a,b,c))
        return a*X + b*Y + c*np.ones(Z.shape) + 0*Z #The +0*Z at the end is to ``re-introduce'' the nans so that the approximation is not plotted everywhere (bit overwhelming)

    def quadratic_appr(self,X,Y,Z):
        Xf,Yf,Zf = X.flatten(), Y.flatten(), Z.flatten()
        no_nans = ~np.isnan(Zf)
        Xf,Yf,Zf = Xf[no_nans],Yf[no_nans],Zf[no_nans]
        A = np.array([Xf*Xf,Xf*Yf,Yf*Yf,Xf,Yf,np.ones(Zf.shape[0])]).T
        a, b, c, d, e, f =  np.linalg.solve(np.matmul(A.T,A),np.matmul(A.T,np.array([Zf]).T)).flatten()
        print('a = %0.4f, b = %0.4f, c = %0.4f, d = %0.4f, e = %0.4f, f = %0.4f'%(a,b,c,d,e,f))
        return a*X*X + b*X*Y + c*Y*Y + d*X + e*Y + f*np.ones(Z.shape) + 0*Z #The +0*Z at the end is to ``re-introduce'' the nans so that the approximation is not plotted everywhere (bit overwhelming)

    def my3Dplot(self):
        H_mesh,C_mesh = np.meshgrid(self.H2,self.C2)

        intrusion_linear = self.linear_appr(H_mesh,C_mesh,self.intrusion_mass)
        intrusion_quadratic = self.quadratic_appr(H_mesh,C_mesh,self.intrusion_mass)

        fig = make_subplots(rows=1, cols=2, 
            specs=[[{'type':'surface'}]*2], 
            subplot_titles=(r'$\text{Intrusion Mass}$', r'$\text{COM }x\text{-coordinate}$'),
            horizontal_spacing=0.15)

        fig.add_trace(go.Surface(z=self.intrusion_mass, x=H_mesh, y=C_mesh, colorscale='viridis',colorbar={'x':0.43,'title':'data'}), row=1, col=1)
        fig.add_trace(go.Surface(z=intrusion_linear, x=H_mesh, y=C_mesh, colorscale='Plotly3', colorbar={'x':0.51, 'title':'Linear approximation'}), row=1, col=1)

        fig.add_trace(go.Surface(z=self.COM_x, x=H_mesh, y=C_mesh, colorscale='plasma'), row=1, col=2)

        fig.update_layout(title_text=r'$\text{Settling speed}: U_s=%0.3f$'%self.U_s,
            height = 600, width = 1200,
            scene ={'xaxis_title':r"h_{2,0}",'yaxis_title':r"c_{2,0}",
                'aspectmode':'manual','aspectratio':dict(x=1, y=1, z=1.),'camera':dict(eye=dict(x=1.1, y=1.98, z=0.66))},
            scene2={'xaxis_title':r'h_{2,0}','yaxis_title':r'c_{2,0}',
                'aspectmode':'manual','aspectratio':dict(x=1, y=1, z=1.),'camera':dict(eye=dict(x=1.1, y=1.98, z=0.66))})
        fig.write_html(self.rootFile + 'solutions/plots/3Dview_' + self.fileName +'.html', include_mathjax="cdn")
        del fig

def make_deposition_plots(US=[0.005,0.01,0.015]):
    '''
    This function runs creates plot when you want to plot multiple surface for differnt settling speeds at once. 
    The DepositionAnalysis class has to be called for each settling speed to generate the data. 
    '''
    def singleSettlingSpeeds():
        for u in US:
            d = DepositionAnalysis(u,'SedimentationInitialConditionTest_2025Jun7/')
            d.my3Dplot()
            d.plot_dimensional_analysis()

    def layeredSettlingSpeeds():
        fig = make_subplots(rows=1, cols=2, 
            specs=[[{'type':'surface'}]*2], 
            subplot_titles=(r'$\text{Intrusion Mass}$', r'$\text{COM }x\text{-coordinate}$'))
  
        settlingSpeedCM = ['Blugrn','Plotly3','Burg']
         
        for i,u in enumerate(US):
            d = DepositionAnalysis(u,'SedimentationInitialConditionTest_2025Jun7/')
            H_mesh,C_mesh = np.meshgrid(d.H2,d.C2)
             

            fig.add_trace(go.Surface(z=d.intrusion_mass, x=H_mesh, y=C_mesh, 
                colorscale=settlingSpeedCM[i],
                colorbar={'x':0.29+0.07*i,'title':{'text':'$U_s=%0.3f$'%(u),'side':'top'},'len':0.75,'thickness':20}), 
                row=1, col=1)
            fig.add_trace(go.Surface(z=d.COM_x, x=H_mesh, y=C_mesh, 
                colorscale=settlingSpeedCM[i],
                colorbar={'x':0.79+0.07*i,'title':{'text':'$U_s=%0.3f$'%(u),'side':'top'},'len':0.75,'thickness':20}), 
                row=1, col=2)

        fig.update_layout(height = 600, width = 1200,
            scene ={'xaxis_title':r'h_{2,0}','yaxis_title':r'c_{2,0}','domain':{'x':[0,.29]},
                'aspectmode':'manual','aspectratio':dict(x=1, y=1, z=1.),'camera':dict(eye=dict(x=1.5, y=2.7, z=0.9))},
            scene2={'xaxis_title':r'h_{2,0}','yaxis_title':r'c_{2,0}','domain':{'x':[0.5,.79]},
                'aspectmode':'manual','aspectratio':dict(x=1, y=1, z=1.),'camera':dict(eye=dict(x=1.5, y=2.7, z=0.9))})
        tempFileName = d.fileName[:d.fileName.find('Us')]+d.fileName[d.fileName.find('Us') + d.fileName[:d.fileName.find('Us')].find('_')+2:]
        fig.write_html(d.rootFile + 'solutions/plots/3Dview_layeredSettling_' + tempFileName + '.html', include_mathjax="cdn")
        del fig
    singleSettlingSpeeds()
    layeredSettlingSpeeds()

def sediment_check(u_s,rootFile,N=5000,sharp=50):
    x = np.linspace(0.7,1.42,73)
    X = np.zeros((x.size,x.size))
    for i in range(x.size):
        for j in range(x.size):
            if x[i]*x[j]<=1:
                X[i,j] = min(LoadSim(x[i],x[j],u_s,rootFile,['d1','d2'],N=N,sharp=sharp).d1.shape[0]-1,1)
    return X

def collision_details(u_s,rootFile,N=5000,sharp=50):
    x = np.linspace(0.7,1.42,73)
    X = np.zeros((x.size,x.size))
    T = np.zeros((x.size,x.size))
    for i in range(x.size):
        for j in range(x.size):
            if x[i]*x[j]<=1:
                X[i,j] = LoadSim(x[i],x[j],u_s,rootFile,[],N=N,sharp=sharp).coll_loc
                T[i,j] = LoadSim(x[i],x[j],u_s,rootFile,[],N=N,sharp=sharp).coll_time
            else:
                X[i,j] = np.nan
                T[i,j] = np.nan
    plt.figure(figsize = [12,5])
    
    plt.rcParams.update({"text.usetex":True})
    figure_title = ['Collision $x$-location','Collision time']
    for i,A in enumerate([X,T]):
        plt.subplot(1,2,i+1)
        pplot = plt.pcolormesh(x,x,A,shading='nearest')
        plt.contour(x,x,A,colors = 'k',linewidths=0.75)
        plt.gca().set_aspect('equal')
        plt.colorbar(pplot)
        plt.title(figure_title[i])
        plt.xlabel('$h_{2,0}$')
        plt.ylabel('$c_{2,0}$')
    plt.tight_layout()
    plt.savefig(rootFile + 'solutions/plots/collision_details_Us%0.3f'%(u_s) + '.pdf')
    plt.show()
    plt.rcParams.update({"text.usetex":False})
