import lc_utils as lcu
from Period_Finding import BLS

qmin = 0.01
qmax = 0.15
output_dir = "/home/echickle/out/"
# data_dir = "/matchfiles/data2/ATLAS/"
data_dir = '/data/ATLAS/WDs/'
wd_main = "/home/echickle/data/GaiaEDR3_WD_main.fits"
rp_ext = "/home/echickle/data/GaiaEDR3_WD_RPM_ext.fits"

gid = 6202621791260175232
fname_atlas = data_dir + str(gid)
suffix="_"+str(gid)

t, y, dy, ra, dec = lcu.load_atlas_lc(fname_atlas, pos_iqr=3, neg_iqr=10)
print(ra)
print(dec)


freqs_to_remove = []
df = 0.05
freqs_to_remove.append([1 - df, 1 + df])
freqs_to_remove.append([1/2. - df, 1/2. + df])
freqs_to_remove.append([1/4. - df, 1/4. + df])


t, y, dy, period, bls_power_best, freqs, power, q, phi0 = \
    BLS(t,y,dy,pmin=2,pmax=10,qmin=qmin,qmax=qmax,freqs_to_remove=freqs_to_remove)
res = lcu.vet_plot(t, y, freqs, power, q, phi0, output_dir=output_dir,
                   objid=gid, objid_type='GAIAID',
                   dy=dy, suffix=suffix, wd_main=wd_main, rp_ext=rp_ext)
sig, snr, wid, period, period_min, q, phi0, epo, rp, nt, dphi = res
print(epo)
print(period*1440)
print(period*q*1440)

# _, _, _, ls_period, ls_freqs, ls_power = LS(t,y,dy,freqs_to_remove=freqs_to_remove)
# suffix='_ra_{}_dec_{}'.format(ra, dec)

# res=lcu.vet_plot(t, y, ls_freqs, ls_power,output_dir=output_dir,
#                  objid=gid, objid_type=None, suffix=suffix, ra=ra,
#                  dec=dec, wd_main=wd_main, rp_ext=rp_ext, 
#                  bls=False)
    
# # per, q, epo = res[3], res[5], res[7]
# # lcu.plot_eclipse_timing(t, y, per, epo, q, output_dir+'GAIAID_{}_{}_{}_'.format(gid, ra, dec))
