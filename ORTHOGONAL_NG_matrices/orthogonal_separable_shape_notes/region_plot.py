import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from scipy import interpolate
from matplotlib.ticker import AutoMinorLocator
from matplotlib.ticker import LogLocator
from matplotlib.ticker import NullFormatter
import warnings

from matplotlib import rc
rc('font',**{'family':'serif','serif':['Computer Modern Roman'],'monospace':['Computer Modern Typewriter']})
rc('text',usetex=True)
mpl.rcParams['text.latex.preamble']=[r'\usepackage{amsmath}',r'\usepackage{amssymb}',r'\usepackage{bm}']
#--> https://stackoverflow.com/questions/23824687/

#@np.vectorize
def S(x1,x2,Delta):
    overall=np.power(9,3*Delta-1)
    A=np.power((x1*x2)/(np.power(1+x1+x2,3)),2*Delta)
    B=np.power(1+np.power(x1,3)+np.power(x2,3),2)
    C=np.power(x1*x2,3)
    out=(overall*A*B)/C
    return out
#--> k3 is equal to one, equivalently we multiply by k3^6...

x1_range=np.linspace(0,1,3000)
x2_range=np.linspace(1,0.5,3000)

@np.vectorize
def bound_plot(x):
    if x<0.5:
        out=1-x
    else:
        out=x
    return out

xpts=x1_range
ypts=bound_plot(xpts)

X,Y=np.meshgrid(x1_range,x2_range)

#**************************************************************************************************************

def col_fun_gl(x):
    return mpl.cm.Blues_r(x)

def col_fun_cfc(x):
    return mpl.cm.Oranges_r(x)

#print("-> plotting...")
#print("\n")

dict_labels=[#"1/2",
"1"]
selected_plot_Delta={"1/2":0.5,"1":1}
#norm={"0":1/3,"1":1/(8/9)}
v_max={"1/2":40,"1":40}
#--> the max one for 0 and 0.5 is completely at random: there I choose the max value of the data, since it diverges...
#v_max={}
#for el in v_min:
#    v_max[el]=???+(???-v_min[el])

save_name={"1/2":"0.5","1":"plot"}
plot_name={"1/2":"1/2","1":"plot"}
#--> kind of like Mathematica's default...

for el in dict_labels:

    print("-> plotting \Delta = "+el+"...")
    print("\n")

    fig,ax=plt.subplots(constrained_layout=True)#(1,1,1),aspect=0.5)
    #ax.set_aspect(0.5)

    #ax.tick_params(direction="in")

    ax.minorticks_on()

    ax.tick_params(which="both",direction="in",pad=6)
    #--> defauls seems to be "major"... Yes, see https://matplotlib.org/api/_as_gen/matplotlib.axes.Axes.tick_params.html...

    plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))

    #xfmt=mpl.ticker.ScalarFormatter()#useMathText=True)
    #xfmt=mpl.ticker.ScalarFormatter(useMathText=True)
    #xfmt.set_powerlimits((-1,1))

    #ax.set_xscale('log')

    xrange=[0,1]
    yrange=[0.5,1]

    ax.set_xlim(xrange[0],xrange[1])
    ax.set_ylim(yrange[0],yrange[1])

    #locmaj = LogLocator(base=10,numticks=15)
    #ax.yaxis.set_major_locator(locmaj)

    #locmin = LogLocator(base=10.0,subs=(0.2,0.4,0.6,0.8),numticks=15)
    #ax.yaxis.set_minor_locator(locmin)
    #ax.yaxis.set_minor_formatter(NullFormatter())

    ax.set_xlabel('$k_1/k_3$',fontsize=19)
    ax.set_ylabel('$k_2/k_3$',fontsize=19)

    #plt.title("$k^3_3(\\langle\\zeta_{{\\bm k}_1}\\zeta_{{\\bm k}_2}\\zeta_{{\\bm k}_3}\\rangle^\\prime)^2/\\prod_{i=1}^3P_\\zeta(k_i)$\\quad ($\\Delta="+plot_name[el]+"$)",y=1.01,fontsize=17) #--> it was 1.025...

    ax.xaxis.set_ticks_position('both')
    ax.yaxis.set_ticks_position('both')

    for tick in ax.xaxis.get_major_ticks():
        tick.label.set_fontsize(17)

    for tick in ax.yaxis.get_major_ticks():
        tick.label.set_fontsize(17)

    #ax.tick_params(direction="in")

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        #Z=norm[el]*S(X,Y,selected_plot_Delta[el])
        Z=S(X,Y,selected_plot_Delta[el])

    #print(Z) #--> for some reason if I had put S = constant this would not work...
    #print(np.shape(x1_range)[0])
    Z=np.full((np.shape(x1_range)[0],np.shape(x2_range)[0]),230)
    Z=np.ma.masked_where(np.logical_or(X+Y<1,X>Y),Z)

    MIN=1#v_min[el]
    MAX=v_max[el]

    print(Z)
    print(np.shape(Z))

    #CS = ax.contourf(Y,X,Z,cmap=plt.cm.plasma,levels=levels[el])#,alpha=0.5)
    #CS=ax.imshow(Z,extent=[X.min(),X.max(),Y.min(),Y.max()],cmap=plt.cm.afmhot)#,vmin=v_min[el])#,vmax=v_max[el])
    CS=ax.imshow(Z,extent=[X.min(),X.max(),Y.min(),Y.max()],cmap='gray',vmin=0,vmax=255)

    #cbar = fig.colorbar(CS,location="bottom",#extend='both',
    #                    shrink=0.9)#,ax=ax) #--> ax should not be needed if later I write cbar.ax.set_xlabel...

    #cbar.set_ticks([1,20,40,60,80,100])
    #cbar.set_ticks([1,10,20,30,40])
    #cbar.ax.tick_params(labelsize=17)

    ax.text(#0.0725
            0.075, 0.96, 'squeezed', style='italic',
                    bbox={'facecolor': 'red', 'alpha': 0.5})#, 'pad': 10})

    ax.text(#0.8125
            0.81, 0.96, 'equilateral', style='italic',
                    bbox={'facecolor': 'red', 'alpha': 0.5})#, 'pad': 10})

    ax.text(0.5, #0.565
            0.57, 'folded', style='italic', horizontalalignment='center',
                    bbox={'facecolor': 'red', 'alpha': 0.5})#, 'pad': 10})

    #cbar.ax.set_xlabel('$k_3^6\\langle\\zeta\\zeta\\zeta\\rangle^\\prime(x_1,x_2,k_3)$')
    #cbar.ax.set_xlabel('$\\langle\\zeta_{x_1}\\zeta_{x_2}\\zeta_{1}\\rangle^\\prime\\,\\,,$\\quad normalized to $\\langle\\zeta_{1}\\zeta_{1}\\zeta_{1}\\rangle^\\prime=1$')
    #cbar.ax.set_xlabel('normalized to $k^3_3(\\langle\\zeta_{{\\bm k}_3}\\zeta_{{\\bm k}_3}\\zeta_{{\\bm k}_3}\\rangle^\\prime)^2/P^3_\\zeta(k_3)=1$',fontsize=17)#,labelpad=1)#,y=0.9)

    ax.plot(xpts,ypts,linestyle="-",color="k",alpha=1,lw=2) #--> can change color and style of this line. Not important for now!!! We will see when/if we need this plot!!!

    plt.savefig('./'+'region_'+save_name[el]+'.pdf', #bbox_extra_artists=(text,),
                bbox_inches='tight')

    plt.close()
