import cuvarbase.bls as bls
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

        tmean = np.mean(t)
        t=t-tmean

        # set up search parameters
        search_params = dict(qmin=qmin, qmax=qmax,
                     # The logarithmic spacing of q
                            # dlogq=0.1,
                            dlogq=0.1,

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
        sols = np.array(sols)
        
        if remove:
                freqs_to_remove = []

                df = 0.1
                freqs_to_remove.append([86400/(200*2) - df, 86400/(200*2) + df])
                freqs_to_remove.append([86400/500 - df, 86400/500 + df])    
                freqs_to_remove.append([86400/(200*3) - df, 86400/(200*3) + df])
                freqs_to_remove.append([86400/600 - df, 86400/600 + df])    
                freqs_to_remove.append([86400/(200*4) - df, 86400/(200*4) + df])
                freqs_to_remove.append([86400/(200*5) - df, 86400/(200*5) + df])     
                freqs_to_remove.append([86400/(200*6) - df, 86400/(200*6) + df]) 
                freqs_to_remove.append([86400/(200*7) - df, 86400/(200*7) + df])   

                for pair in freqs_to_remove:
                        idx = np.where((freqs < pair[0]) | (freqs > pair[1]))[0]
                        freqs = freqs[idx]
                        bls_power = bls_power[idx]
                        sols = sols[idx]

        i_best = np.argmax(bls_power)
        bls_power_best = bls_power[i_best]
        f_best = freqs[i_best]
        period=1.0/f_best
        q_best, phi0_best = sols[i_best]   

        # >> finding second peak        
        # freqs_to_remove = [[f_best-0.1*f_best, f_best+0.1*f_best]]
        # for pair in freqs_to_remove:
        #         idx = np.where((freqs < pair[0]) | (freqs > pair[1]))[0]
        #         bls_power2 = bls_power[idx]
        #         freqs2 = freqs[idx]
        # bls_power_best=(np.max(bls_power)-np.mean(bls_power2))/np.std(bls_power2)        
        t=t+tmean # >> add back the mean 

        return t, y, dy, period, bls_power_best, freqs, bls_power, q_best, phi0_best

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

def LS_Astropy(t,y,dy, remove=True,pmax=0.25):
        from astropy.timeseries import LombScargle
        t=t-np.mean(t)
        # lightcurves=[(t,y,dy)]
        # proc = gls.LombScargleAsyncProcess()
        # freqs=np.linspace(1/0.25,86400/400,100000) 
        # result = proc.run([(t, y, dy)], freqs=freqs)
        # proc.finish()
        # freqs, ls_power = result[0]
        freqs, ls_power = LombScargle(t,y).autopower()
        
        if remove:
                freqs_to_remove = []
                df = 0.1
                freqs_to_remove.append([86400/(200*2) - df, 86400/(200*2) + df])
                freqs_to_remove.append([86400/500 - df, 86400/500 + df])    
                freqs_to_remove.append([86400/(200*3) - df, 86400/(200*3) + df])
                freqs_to_remove.append([86400/600 - df, 86400/600 + df])    
                freqs_to_remove.append([86400/(200*4) - df, 86400/(200*4) + df])
                freqs_to_remove.append([86400/(200*5) - df, 86400/(200*5) + df])     
                freqs_to_remove.append([86400/(200*6) - df, 86400/(200*6) + df]) 
                freqs_to_remove.append([86400/(200*7) - df, 86400/(200*7) + df])   
                for pair in freqs_to_remove:
                        idx = np.where((freqs < pair[0]) | (freqs > pair[1]))[0]
                        freqs = freqs[idx]
                        ls_power = ls_power[idx]
                idx = np.where(freqs < 86400/400)[0]
                freqs = freqs[idx]
                ls_power = ls_power[idx]
                idx = np.where(freqs > 1/pmax)[0]
                freqs = freqs[idx]
                ls_power = ls_power[idx]                                        
        best_freq=freqs[np.argmax(ls_power)]
        significance=(np.max(ls_power)-np.mean(ls_power))/np.std(ls_power)      
        period=1.0/freqs[np.argmax(ls_power)]
        return t, y, dy, period, significance, freqs, ls_power        

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

