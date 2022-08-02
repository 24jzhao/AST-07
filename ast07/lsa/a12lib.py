#required vals to be defined:
#oname, av1, obj, an1 (obj = an1), object_number, some other stuff
import matplotlib.pyplot as plt
from itertools import chain
from matplotlib.backends.backend_pdf import PdfPages
def data_format (object_number, all_mags, all_mjd, all_magerrs):
    filts = []
    temp_u = ['u' for i in all_mags[object_number][0]]
    temp_g = ['g' for i in all_mags[object_number][1]]
    temp_r = ['r' for i in all_mags[object_number][2]]
    temp_i = ['i' for i in all_mags[object_number][3]]
    temp_z = ['z' for i in all_mags[object_number][4]]
    filts.append(temp_u)
    filts.append(temp_g)
    filts.append(temp_r)
    filts.append(temp_i)
    filts.append(temp_z)
    filts = list(chain.from_iterable(filts))
    mags = list(chain.from_iterable(all_mags[object_number])) #turning 2d array into 1d
    t = list(chain.from_iterable(all_mjd[object_number]))
    #dy = [1 for i in mags]
    dy = list(chain.from_iterable(all_magerrs[object_number]))
    dy = [99 if i == 0 else i for i in dy]
    return t, mags, dy, filts


def folded_light_curve(obj, period, pdf_name, oname, av1, all_mjd, all_mags):
    plt.ioff()
    with PdfPages(pdf_name) as pdf:
        plt.figure(figsize = (9, 12))
        plt.xlabel('Phase')
        plt.ylabel('Magnitude')
        plt.gca().invert_yaxis()
        u_phase = [(i%period)/period for i in all_mjd[obj][0]]
        g_phase = [(i%period)/period for i in all_mjd[obj][1]]
        r_phase = [(i%period)/period for i in all_mjd[obj][2]]
        i_phase = [(i%period)/period for i in all_mjd[obj][3]]
        z_phase = [(i%period)/period for i in all_mjd[obj][4]]
        #f*p^2/length of observations, f is about 0.1, f = phase error for individual cycle
        plt.scatter(u_phase, all_mags[obj][0], s = 5, c = 'blue', label = 'u')
        plt.scatter(g_phase, all_mags[obj][1], s = 5, c = 'green', label = 'g')
        plt.scatter(r_phase, all_mags[obj][2], s = 5, c = 'purple', label = 'r')
        plt.scatter(i_phase, all_mags[obj][3], s = 5, c = 'gold', label = 'i')
        plt.scatter(z_phase, all_mags[obj][4], s = 5, c = 'tab:red', label = 'z')
        plt.legend()
        
        #oname is a dynamic 1d list containing global index and object name (in that order).
        plt.title(str(oname[1])+ " (" + str(av1) + ")") 
        
        #TODO: map new object names without destroying it 
        pdf.savefig()
        plt.close()

