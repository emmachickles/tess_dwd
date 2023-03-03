import cuvarbase.bls as bls
import cuvarbase.lombscargle as gls
import numpy as np

def frequency_grid(t,y,pmin=3,pmax=True,qmin=2e-2,qmax=0.12):
        baseline = max(t) - min(t)
        df = qmin / baseline
        if type(pmax) == type(True):
                fmin = 4/baseline
        else:
                fmin=1/pmax
        fmax = 1440.0/pmin # pmin in minutes
        nf = int(np.ceil((fmax - fmin) / df))
        freqs = fmin + df * np.arange(nf)
        return freqs

def remove_harmonics(freqs, power, dur=None, epo=None):
        freqs_to_remove = []
        for i in [2,3,4,5,6]:
                centr = 86400 / (200*i)
                freqs_to_remove.append([centr - 0.0001, centr + 0.0001])
        for pair in freqs_to_remove:
                idx = np.where((freqs < pair[0]) | (freqs > pair[1]))[0]
                freqs = freqs[idx]
                power = power[idx]
                if dur is not None:
                        dur = dur[idx]
                        epo = epo[idx]
        if dur is not None:
                return freqs, power, dur, epo
        else: 
                return freqs, power
        
def BLS(t,y,dy,pmin=3,pmax=True,qmin=2e-2,qmax=0.12,remove=True):

        t=t-np.mean(t)
        y=y
        dy=dy

        # set up search parameters
        search_params = dict(qmin=qmin, qmax=qmax,
                     # The logarithmic spacing of q
                            # dlogq=0.1,
                            dlogq=0.5,

                     # Number of overlapping phase bins
                     # to use for finding the best phi0
                        noverlap=3)
        
        # derive baseline from the data for consistency
        baseline = max(t) - min(t)

        # df ~ qmin / baseline
        df = search_params['qmin'] / baseline
        if type(pmax) == type(True):
                fmin = 4/baseline
        else:
                fmin=1/pmax

        fmax = 1440.0/pmin # pmin in minutes

        nf = int(np.ceil((fmax - fmin) / df))
        freqs = fmin + df * np.arange(nf)
        
        # bls_power = bls.eebls_gpu_fast(t, y, dy, freqs,
        #                                **search_params)
        bls_power, sols = bls.eebls_gpu(t, y, dy, freqs,
                                       **search_params)

        
        # if remove:
        #         # low=142
        #         # high=146
        #         # freqs_to_remove = [[low,high], [low/2,high/2],[low/3,high/3], [low/4,high/4], [low/5,high/5], [low/6,high/6], [low/7,high/7],[low/8,high/8], [low/9,high/9]]

        #         rmv = 10.
        #         dp = 0.4
        #         pers_to_remove = []
        #         for i in range(1,10):
        #                 # pers_to_remove.append([rmv*i - dp*i, rmv*i + dp*i])
        #                 pers_to_remove.append([rmv*i - dp, rmv*i + dp])

        #         idx = np.empty(0, dtype='int')
        #         for pair in pers_to_remove:
        #                 idx = np.append(idx, np.where((freqs < 1440/pair[1]) | (freqs > 1440/pair[0]))[0])

        #         f_best = freqs[idx][np.argmax(bls_power[idx])]
        #         p_best = np.max(bls_power[idx])

        #         # >> check if there is a peak at multiples
        #         f1 = f_best/2
        #         f2 = f_best*2
        #         for f in [f1, f2]:
        #                 if f > np.min(freqs) and f < np.max(freqs):
        #                         new_ind = np.argsort( np.abs( freqs - f1  ) )[0]
        #                         new_pow = bls_power[new_ind]
        #                         if new_pow > p_best:
        #                                 f_best = f
        #                                 p_best = new_pow
        #         # import pdb
        #         # pdb.set_trace()
        # if remove:
        #         freq=360
        #         delta=0.025
        #         freqs_to_remove=[]
        #         n=1
        #         while n>0:
        #                 freqs_to_remove.append([freq*n-delta/(n/5)**0.5,freq*n+delta/(n/5)**0.5])
        #                 n=n-0.02
        #                 #freqs_to_remove = [[3e-2,4e-2], [49.99,50.01], [48.99,49.01], [47.99,48.01], [46.99,47.01], [45.99,46.01], [3.95,4.05], [2.95,3.05], [1.95,2.05], [0.95,1.05], [0.48, 0.52], [0.32, 0.34], [0.24, 0.26], [0.19, 0.21]]
        #         for pair in freqs_to_remove:
        #                 idx = np.where((freqs < pair[0]) | (freqs > pair[1]))[0]
        #                 freqs = freqs[idx]
        #                 bls_power = bls_power[idx]


        # fname = open('/home/submit/echickle/foo.txt', 'a')
        
        if remove:
                freqs_to_remove = []
                # freqs_to_remove.append([86400/400 - 0.1, 86400/400 + 0.1])
                # freqs_to_remove.append([86400/600 - 0.07, 86400/600 + 0.07])
                # freqs_to_remove.append([86400/800 - 0.03, 86400/800 + 0.03])
 
                freqs_to_remove.append([86400/(200*2) - 1.2, 86400/(200*2) + 1.2])
                freqs_to_remove.append([86400/500 - 1, 86400/500 + 1])                
                freqs_to_remove.append([86400/(200*3) - 0.1, 86400/(200*3) + 0.1])
                freqs_to_remove.append([86400/600 - 1, 86400/600 + 1])    
                freqs_to_remove.append([86400/(200*4) - 0.1, 86400/(200*4) + 0.1])
                freqs_to_remove.append([86400/(200*5) - 3, 86400/(200*5) + 3])     
                freqs_to_remove.append([86400/(200*6) - 3, 86400/(200*6) + 3]) 
                freqs_to_remove.append([86400/(200*7) - 2, 86400/(200*7) + 2])               
                # for i in [2,3,4]:
                #         centr = 86400 / (200*i)
                #         freqs_to_remove.append([centr - 0.05, centr + 0.05])
                        # fname.write(str(i)+' '+str(centr)+'\n')
                        # fname.write(str(i)+' '+str([centr - 0.0001, centr + 0.0001])+'\n')
                for pair in freqs_to_remove:
                        idx = np.where((freqs < pair[0]) | (freqs > pair[1]))[0]
                        freqs = freqs[idx]
                        bls_power = bls_power[idx]
        
        f_best = freqs[np.argmax(bls_power)]
        # fname.write('f_best '+str(f_best)+'\n\n')
        # fname.close()
        period=1.0/f_best
        q_best, phi0_best = sols[np.argmax(bls_power)]        

        # -- significance of peak ----------------------------------------------
        # >> compare to median
        bls_power_best=(np.max(bls_power)-np.median(bls_power))/(np.std(bls_power))
        
        # >> finding second peak        
        # freqs_to_remove = [[f_best-0.1*f_best, f_best+0.1*f_best]]
        # for pair in freqs_to_remove:
        #         idx = np.where((freqs < pair[0]) | (freqs > pair[1]))[0]
        #         bls_power2 = bls_power[idx]
        #         freqs2 = freqs[idx]
        # bls_power_best=(np.max(bls_power)-np.mean(bls_power2))/np.std(bls_power2)        

        return t, y, dy, period, bls_power_best, freqs, bls_power, q_best, phi0_best

# >> seems to be a repeat of the above function
# def BLS_Full(t,y,dy,pmin=3,pmax=True,qmin=2e-2,qmax=0.12,remove=True,trim=True):

#         t=t-np.mean(t)
#         y=y
#         dy=dy

# # generate data with a transit


# # set up search parameters
#         search_params = dict(qmin=qmin, qmax=qmax,
#                      # The logarithmic spacing of q
#                             dlogq=0.1,

#                      # Number of overlapping phase bins
#                      # to use for finding the best phi0
#                         noverlap=3)
        
# # derive baseline from the data for consistency
#         baseline = max(t) - min(t)

# # df ~ qmin / baseline
#         df = search_params['qmin'] / baseline
#         if pmax:
#                 fmin = 4/baseline
#         else:
#                 fmin=1/pmax

#         fmax = 1440.0/pmin

        
#         nf = int(np.ceil((fmax - fmin) / df))
#         freqs = fmin + df * np.arange(nf)
#         if remove:
#                 low=142
#                 high=146
#                 freqs_to_remove = [[low,high], [low/2,high/2],[low/3,high/3], [low/4,high/4], [low/5,high/5], [low/6,high/6], [low/7,high/7],[low/8,high/8], [low/9,high/9]]
#                 #freqs_to_remove = [[3e-2,4e-2], [49.99,50.01], [48.99,49.01], [47.99,48.01], [46.99,47.01], [45.99,46.01], [3.95,4.05], [2.95,3.05], [1.95,2.05], [0.95,1.05], [0.48, 0.52], [0.32, 0.34], [0.24, 0.26], [0.19, 0.21]]
#                 for pair in freqs_to_remove:
#                         idx = np.where((freqs < pair[0]) | (freqs > pair[1]))[0]
#                         freqs = freqs[idx]

#         bls_power = bls.eebls_gpu_fast(t, y, dy, freqs,
#                                    **search_params)
#         f_best = freqs[np.argmax(bls_power)]                               
#         if trim:       

#                 best_freq=freqs[np.argmax(bls_power)]
#                 freqs_to_remove = [[best_freq-0.1*best_freq, best_freq+0.1*best_freq]]
#                 for pair in freqs_to_remove:
#                         idx = np.where((freqs < pair[0]) | (freqs > pair[1]))[0]
#                         bls_power2 = bls_power[idx]
#                         freqs2 = freqs[idx]
                        


                
                                        
#                 bls_power_best=(np.max(bls_power)-np.median(bls_power2))/np.std(bls_power2)
#         else:


#                 bls_power_best=(np.max(bls_power)-np.median(bls_power))/(np.std(bls_power))
#         period=1.0/f_best


#         return t, y, dy, period, bls_power_best, freqs, bls_power

def LS(t,y,dy,pmin=2,remove=False):
        t=t-np.mean(t)
        lightcurves=[(t,y,dy)]
        baseline=np.max(t)-np.min(t)

        df = 1.0 / (baseline*3.0)
        fmin = 4.0/baseline
        fmax = 1440/pmin

        nf = int(np.ceil((fmax - fmin) / df))
        freqs=np.linspace(fmin,fmax,nf)

        proc = gls.LombScargleAsyncProcess(use_double=True)
        result = proc.run([(t, y, dy)],freqs=freqs,use_fft=True)
        proc.finish()
        freqs, ls_power = result[0]

        if remove:
                freqs_to_remove = [[3e-2,4e-2], [49.99,50.01], [48.99,49.01], [47.99,48.01], [46.99,47.01], [45.99,46.01], [3.95,4.05], [2.95,3.05], [1.95,2.05], [0.95,1.05], [0.48, 0.52], [0.32, 0.34], [0.24, 0.26], [0.19, 0.21]]
                for pair in freqs_to_remove:
                        idx = np.where((freqs < pair[0]) | (freqs > pair[1]))[0]
                        ls_power = ls_power[idx]
                        freqs = freqs[idx]
                                
        significance=(np.max(ls_power)-np.mean(ls_power))/np.std(ls_power)

        period=1.0/freqs[np.argmax(ls_power)]

        return t, y, dy, period, significance
        
def LS_Full(t,y,dy,pmin=2,pmax=True,oversample_factor=3.0,trim=False):
        t=t-np.mean(t)
        lightcurves=[(t,y,dy)]
        baseline=np.max(t)-np.min(t)

        df = 1.0 / (baseline*oversample_factor)
        if type(pmax) == type(True):
                fmin = 4/baseline
        else:
                fmin=1/pmax                
        fmax = 1440/pmin

        nf = int(np.ceil((fmax - fmin) / df))
        freqs=np.linspace(fmin,fmax,nf)

        proc = gls.LombScargleAsyncProcess(use_double=True)
        result = proc.run([(t, y, dy)],freqs=freqs,use_fft=True)
        proc.finish()
        freqs, ls_power = result[0]
        #freqs_to_remove = [[3e-2,4e-2], [49.99,50.01], [48.99,49.01], [47.99,48.01], [46.99,47.01], [45.99,46.01], [3.95,4.05], [2.95,3.05], [1.95,2.05], [0.95,1.05], [0.48, 0.52], [0.32, 0.34], [0.24, 0.26], [0.19, 0.21]]
        #for pair in freqs_to_remove:
        #       idx = np.where((freqs < pair[0]) | (freqs > pair[1]))[0]
        #       ls_power = ls_power[idx]
        #       freqs = freqs[idx]
        
        best_freq=freqs[np.argmax(ls_power)]
        if trim:
                freqs_to_remove = [[best_freq-0.1*best_freq, best_freq+0.1*best_freq]]
                for pair in freqs_to_remove:
                        idx = np.where((freqs < pair[0]) | (freqs > pair[1]))[0]
                        ls_power2 = ls_power[idx]
                        freqs2 = freqs[idx]
                        

                f_best = freqs[np.argmax(ls_power)]
                
                                        
                significance=(np.max(ls_power)-np.mean(ls_power2))/np.std(ls_power2)
        else:
        
                f_best = freqs[np.argmax(ls_power)]
                
                                        
                significance=(np.max(ls_power)-np.mean(ls_power))/np.std(ls_power)      

        period=1.0/freqs[np.argmax(ls_power)]

        return t, y, dy, period, significance, freqs, ls_power

