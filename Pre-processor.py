import numpy as np
import matplotlib.pyplot as plt
from os import listdir
import imageio
import os

#%%

# Data paths (delete as appropriate)
# MVSC (Windows users):
path = './x64/Release/'
# Linux/iOS users:
#path = './Life/'

# ONLY NECESSARY IF YOU WISH TO USE AN IMAGE
use_image = True        # If true, please fill in the image_path
image_path = './Images/starry_night.png' 
# Set the number of cores you will run the Game of Life on
p = 14    # It is vital that this matches the runtime command (e.g. `mpiexec -n p Life.exe`)
# Set the brightness ratio (higher => darker)
threshold = 0.4
domain_size = 100


#%%

# Firstly, clean up all directories that will be written into throught the pipeline of the program run
# Delete files from folder that contains input binary digits to fill in the Game of Life local grids
filelist = [ f for f in os.listdir(path+'indata/') if f.endswith(".txt") ]
for f in filelist:
    os.remove(os.path.join(path+'indata/', f))
# Delete files from folder that contains meta data to be used by post-processor
filelist = [ f for f in os.listdir(path+'meta/') if f.endswith(".txt") ]
for f in filelist:
    os.remove(os.path.join(path+'meta/', f))
# Delete files from folder that contains output grid data to be used by post-processor
filelist = [ f for f in os.listdir(path+'data/') if f.endswith(".txt") ]
for f in filelist:
    os.remove(os.path.join(path+'data/', f))

if (use_image == True):
    # Create a dense image to play the Game of Life upon
    # Load image as grayscale
    dense_image = np.array(imageio.imread(image_path, as_gray=True) / 255)
    # We want it to be roughly (100 x 100)
    len_i, len_j = np.shape(dense_image)
    scale_i = int(len_i / 100)
    scale_j = int(len_j / 100)

    # Downscale it
    image = dense_image[::scale_i, ::scale_j]
    rows_global, cols_global = np.shape(image)

    # Do domain decomposition to calculate how to partition the domain
    if (p > 1):
        gap = p
        for n in np.arange(1,p):
            if (p % n == 0):
                if (gap > abs(n - p/n)):
                    gap = abs(n - p/n);
                    proc_i = n
                    proc_j = int(p / n)
    else:
        proc_i = 1
        proc_j = 1

    # Preallocate memory to store the id and how many rows and columns each processor is responsible for
    id_array = np.arange(p, dtype=int)
    imax_array = np.zeros(p)
    jmax_array = np.zeros(p)

    # Calculate position of processor in domain
    i_id = (np.floor(id_array / proc_j)).astype(int)
    j_id = (id_array % proc_j).astype(int)

    # Calculate how to distribute rows amongst processors
    for me in range(p):
        rows_rem = rows_global
        for i in range(i_id[me] + 1):
            imax_array[me] = int(rows_rem / (proc_i - i))
            rows_rem -= imax_array[me]
    # Calculate how to distribute columns amongst processors 
    for me in range(p):
        cols_rem = cols_global
        for j in range(j_id[me] + 1):
            jmax_array[me] = int(cols_rem / (proc_j - j))
            cols_rem -= jmax_array[me]

    # Partition the image into ones and zeros
    image[image > threshold] = 1
    image[image <= threshold] = 0
    # Convert to booleans
    image.astype(bool)
    row_start = int(0)
    cnt = 0
    for i in range(proc_i):
        row_end = int(np.cumsum(imax_array[::proc_i])[i])
        col_start = int(0)
        for j in range(proc_j):
            col_end = int(np.cumsum(jmax_array[0:proc_j])[j])
            grid = image[row_start:row_end, col_start:col_end].astype(bool)    
            plt.plot([col_start, col_end], [row_end, row_end], 'b-', linewidth = 1)
            plt.plot([col_end, col_end], [row_start, row_end], 'b-', linewidth = 1)

            np.savetxt(path+'indata/grid_%i_%i_%i.txt' % (id_array[cnt], i, j), [grid.ravel()], delimiter = ",", fmt="%d")
            col_start = col_end
            cnt += 1
        row_start = row_end

    plt.imshow(image, 'binary')
    plt.title("Starry Night - Van Gogh")
    plt.show()

    print("Please enter these values into the global variables of the Life.cpp script")
    print("(imax_global, jmax_global) = (%i, %i)" % np.shape(image))

