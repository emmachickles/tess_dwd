import matplotlib.pyplot as plt
import numpy as np
import os


snr_threshold = 10

in_dir  = "/scratch/echickle/s0063/s0063-bls-2-1-230403/"
out_dir = "/scratch/echickle/s0063/s0063-bls-2-1-hsnr/"
out_fig = "/scratch/echickle/s0063/s0063-bls-2-1-snr"
os.makedirs(out_dir)

fnames = os.listdir(in_dir)
snr_list = []
for i in range(len(fnames)):
    snr = float(fnames[i].split('_')[3])
    snr_list.append(snr)
    if snr > snr_threshold and fnames[i].split('.')[-1] == 'png' \
       and snr < 150:
        os.system('cp '+in_dir+fnames[i]+' '+out_dir+fnames[i])
snr_list = np.array(snr_list)
plt.figure()
_ = plt.hist(snr_list, bins=100)
plt.xlabel("SNR")
plt.savefig(out_fig+".png")
print(out_fig+".png")

inds = np.nonzero(snr_list < 100)
plt.figure()
_ = plt.hist(snr_list[inds], bins=100)
plt.xlabel("SNR")
plt.savefig(out_fig+"_zoom.png")
print(out_fig+"_zoom.png")
