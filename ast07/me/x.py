
t, mags, dy, filts = data_format(26059)

#import numpy as np
#import matplotlib.pyplot as plt
#import matplotlib as mpl
#mpl.style.use('ggplot')
#mpl.rc('axes', color_cycle=["#4C72B0", "#55A868", "#C44E52","#8172B2", "#CCB974"]) colors for each band
#from gatspy.periodic import LombScargleMultiband
#from gatspy import datasets, periodic

# Choose a Sesar 2010 object to base our fits on
#lcid = 1019544
#rrlyrae = datasets.RRLyraeGenerated(lcid, random_state=0)

# Generate data in a 6-month observing season
#Nobs = 60
#rng = np.random.RandomState(0)

#nights = np.arange(180)
#rng.shuffle(nights)
#nights = nights[:Nobs]

#t = 57000 + nights + 0.05 * rng.randn(Nobs)
#dy = 0.06 + 0.01 * rng.randn(Nobs)
#mags = np.array([rrlyrae.generated(band, t, err=dy, corrected=False)for band in 'ugriz'])
"""
All the above are not necessary for the NGVS_legacy work, becuase you do have realistic data.
I generate mock datesets just to test to code and show an example result.
"""

"""
Here start the fitting process. 
t, mags, dy, filts are all N-d arraies of your observation data. They mean observation time in MJD, apperent magnitudes, errors, and filter list, respectively
If you do not have error list, set dy to be all ones.
filts could be anything in the format of np.array(['u','u','g', ..., 'z']), as long as its elements are filter names and its length equal to the mags array
"""
#example of appropriate format
#filts = np.take(list('ugriz'), np.arange(Nobs), mode='wrap') 
# 
#mags = mags[np.arange(Nobs) % 5, np.arange(Nobs)]
#masks = [(filts == band) for band in 'ugriz']#separate ugriz to 5 sublists

#below is the necessary code
periods = np.linspace(0.1, 1.0, 100000) # This defines the search range of your period, you can specify it at your will. These are in days.


#2 different ways of fitting
##model = periodic.NaiveMultiband(BaseModel=periodic.LombScargleFast) 
# specify the method to be the naive multiband LS, which means you fit data in each band separately, and get a score list for each band.
# serves as a good first try on your NGVS_legacy data
#above is the fastest way. x axis is period, y axis is score. 5d array output. real variable should have same peak for all bands
##model.fit(t, mags, dy, filts) 
##P = model.scores(periods) 
# This is the fitting score list you want. 
# It is a 5xN array, where N is number of periods tried, P[0] is the fit socres of periods with your u band data. And so on for P[1] for g, P[2] for i, ...
# Each element in P[i] correspond to a period in the array periods you input. The closer to 1, the better.


#all bands done together
LS_multi = periodic.LombScargleMultiband(Nterms_base=1, Nterms_band=0)#initiate structure variable
LS_multi.fit(t, mags, dy, filts)#input our data
P_multi = LS_multi.periodogram(periods)#function where input is periods
# A non-naive way of multiband fitting. This time all data from all bands are fitted simultaneously, means you do not get scores for each band separately.
# P_multi will be a 1-d array, has equal length to your input periods. The maximum value in P_multi will be the best fit, and its corresponding period will be the best period.

#do both!!

"""
From here are visualization of the results.

fig = plt.figure(figsize=(10, 4))
gs = plt.GridSpec(5, 2, left=0.07, right=0.95, bottom=0.15,
                  wspace=0.1, hspace=0.6)
ax = [fig.add_subplot(gs[:, 0]),
      fig.add_subplot(gs[:-2, 1]),
      fig.add_subplot(gs[-2:, 1])]

for band, mask in zip('ugriz', masks):
    ax[0].errorbar((t[mask] / rrlyrae.period) % 1, mags[mask], dy[mask],
                   fmt='.', label=band)
ax[0].set_ylim(18, 14.5)
ax[0].legend(loc='upper left', fontsize=12, ncol=3)
ax[0].set_title('Folded Data, 1 band per night (P={0:.3f} days)'
                ''.format(rrlyrae.period), fontsize=12)
ax[0].set_xlabel('phase')
ax[0].set_ylabel('magnitude')

for i, band in enumerate('ugriz'):
    offset = 4 - i
    ax[1].plot(periods, P[band] + offset, lw=1)
    ax[1].text(0.89, 1 + offset, band, fontsize=10, ha='right', va='top')
ax[1].set_title('Standard Periodogram in Each Band', fontsize=12)
ax[1].yaxis.set_major_formatter(plt.NullFormatter())
ax[1].xaxis.set_major_formatter(plt.NullFormatter())
ax[1].set_ylabel('power + offset')


ax[2].plot(periods, P_multi, lw=1, color='gray')

ax[2].set_title('Multiband Periodogram', fontsize=12)
ax[2].set_yticks([0, 0.5, 1.0])
ax[2].set_ylim(0, 1.0)
ax[2].yaxis.set_major_formatter(plt.NullFormatter())
ax[2].set_xlabel('Period (days)')
ax[2].set_ylabel('power')
"""


plt.ion() #turns figure display on
plt.figure()
plt.scatter(periods, P_multi, s = 0.05)

# |%%--%%| <OEulKEDgbl|HfMQ9k54UZ>

best_period = max(P_multi)
for i in range(len(P_multi)):
    if P_multi[i] == best_period:
        index = i
print(periods[index])

# |%%--%%| <HfMQ9k54UZ|EL5F4xxt1S>

def not_pdf_folded_light_curve(obj, period):
    #rr lyrae has variation in g band of 0.6 or 0.5 in unfolded light curve
    plt.ion()
    #with PdfPages(pdf_name) as pdf:
        #for star_object in objects:
    plt.figure(figsize = (9, 12))
    plt.xlabel('Phase')
    plt.ylabel('Magnitude')
    plt.gca().invert_yaxis()
    peak_mag = min(all_mags[obj][0])
    u_phase = [(i%period)/period for i in all_mjd[obj][0]]
    g_phase = [(i%period)/period for i in all_mjd[obj][1]]
    r_phase = [(i%period)/period for i in all_mjd[obj][2]]
    i_phase = [(i%period)/period for i in all_mjd[obj][3]]
    z_phase = [(i%period)/period for i in all_mjd[obj][4]]
    for i in range(len(all_mags[obj][0])):
        if all_mags[obj][0][i] == peak_mag:
            peak_index = i
    u_phase_change = [u_phase[i]- u_phase[peak_index] if u_phase[i]- u_phase[peak_index] >= 0 else u_phase[i]- u_phase[peak_index] + 1 for i in range(len(u_phase))]
    g_phase_change = [g_phase[i]- u_phase[peak_index] if g_phase[i]- u_phase[peak_index] >= 0 else g_phase[i]- u_phase[peak_index] + 1 for i in range(len(g_phase))]
    r_phase_change = [r_phase[i]- u_phase[peak_index] if r_phase[i]- u_phase[peak_index] >= 0 else r_phase[i]- u_phase[peak_index] + 1 for i in range(len(r_phase))]
    i_phase_change = [i_phase[i]- u_phase[peak_index] if i_phase[i]- u_phase[peak_index] >= 0 else i_phase[i]- u_phase[peak_index] + 1 for i in range(len(i_phase))]
    z_phase_change = [z_phase[i]- u_phase[peak_index] if z_phase[i]- u_phase[peak_index] >= 0 else z_phase[i]- u_phase[peak_index] + 1 for i in range(len(z_phase))]
        #f*p^2/length of observations, f is about 0.1, f = phase error for individual cycle
    plt.scatter(u_phase_change, all_mags[obj][0], s = 5, c = 'blue', label = 'u')
    plt.scatter(g_phase_change, all_mags[obj][1], s = 5, c = 'green', label = 'g')
    plt.scatter(r_phase_change, all_mags[obj][2], s = 5, c = 'purple', label = 'r')
    plt.scatter(i_phase_change, all_mags[obj][3], s = 5, c = 'gold', label = 'i')
    plt.scatter(z_phase_change, all_mags[obj][4], s = 5, c = 'tab:red', label = 'z')
    plt.legend()
    plt.title(fnames[obj] + " (" + str(obj) + ")")
        #pdf.savefig()
        #plt.close()
    #0.6 seconds ideal step
    # i/p_true + n where n is an integer all reciprocated is the beat frequency

# |%%--%%| <EL5F4xxt1S|7uLgHexxQp>

not_pdf_folded_light_curve(26059, periods[index])
