import matplotlib.pyplot as plt
import numpy as np
import os


mse_threshold = 0.02

# in_dir  = "/scratch/echickle/s0063/s0063-ls-3-1-230403/"
# out_dir = "/scratch/echickle/s0063/s0063-ls-3-1-hmse/"
# out_fig = "/scratch/echickle/s0063/s0063-ls-3-1-mse"
in_dir  = "/home/echickle/out/s0057/s0057-ls-2-3-230403/"
out_dir = "/home/echickle/out/s0057/s0057-ls-2-3-lmse/"
out_fig = "/home/echickle/out/s0057/s0057-ls-2-3-mse"

os.makedirs(out_dir)

fnames = os.listdir(in_dir)
mse_list = []
for i in range(len(fnames)):
    mse = float(fnames[i].split('_')[3])
    mse_list.append(mse)
    if mse < mse_threshold and fnames[i].split('.')[-1] == 'png':
        os.system('cp '+in_dir+fnames[i]+' '+out_dir+fnames[i])
mse_list = np.array(mse_list)
plt.figure()
_ = plt.hist(mse_list, bins=100)
plt.xlabel("Mean Squared Error of sine fit")
plt.savefig(out_fig+".png")
print(out_fig+".png")

inds = np.nonzero(mse_list < 1)
plt.figure()
_ = plt.hist(mse_list[inds], bins=100)
plt.xlabel("Mean Squared Error of sine fit")
plt.xscale('log')
plt.savefig(out_fig+"_zoom.png")
print(out_fig+"_zoom.png")
