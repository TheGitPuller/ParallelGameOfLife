import numpy as np
import matplotlib.pyplot as plt
from os import listdir
import matplotlib.animation as anim

#%%

# Data paths (delete as appropriate)
# MVSC (Windows users):
path = './x64/Release/'
# Linux/iOS users:
#path = './Life/'

#%%

# Initialise the number of processors to occupy domain in each direction
n_proc_rows = 0
n_proc_cols = 0
# Extract macro-info using file titles
for file in listdir(path + 'meta/'):
    i = int(file.split('_')[0])
    j = int(file.split('_')[1])
    if (i > n_proc_rows):
        n_proc_rows = i;
    if (j > n_proc_cols):
        n_proc_cols = j;
# Make number inclusive        
n_proc_rows += 1;
n_proc_cols += 1;
n_procs = n_proc_rows*n_proc_cols
# Infer number of iterations from metadata
iter = int(listdir(path + 'meta/')[0].split('_')[2])
# Initialise arrays to signify what processor sits where in the domain
i_array = np.zeros((n_proc_rows * n_proc_cols), dtype=int)
j_array = np.zeros((n_proc_rows * n_proc_cols), dtype=int)
# Initialise arrays to signify 
imax_local_array = np.zeros((n_proc_rows * n_proc_cols), dtype=int)
jmax_local_array = np.zeros((n_proc_rows * n_proc_cols), dtype=int)
time_array = np.zeros((n_proc_rows * n_proc_cols))

for r in range(n_proc_rows):
    for c in range(n_proc_cols):
        string = np.genfromtxt(path + 'meta/%i_%i_%i_%i_info.txt' % (r,c,iter, n_procs))
        i_array[r * n_proc_cols + c] = r
        j_array[r * n_proc_cols + c] = c
        imax_local_array[r * n_proc_cols + c] = int(string[0])
        jmax_local_array[r * n_proc_cols + c] = int(string[1])
        time_array[r * n_proc_cols + c] = float(string[2])

row_cumsum = np.cumsum(imax_local_array[::n_proc_cols])
col_cumsum = np.cumsum(jmax_local_array[0:n_proc_cols])

imax_global = int(row_cumsum[-1])
jmax_global = int(col_cumsum[-1])
matrix = np.zeros((imax_global , jmax_global))
matrix_store = np.zeros((iter, imax_global, jmax_global))

for it in range(iter):
    row_start = 0
    line_artist = []
    for r in range(n_proc_rows):
        row_end = row_cumsum[r]
        col_start = 0
        for c in range(n_proc_cols):
            A = np.genfromtxt(path + 'data/%i_%i_%i.txt' % (r,c,it))
            index = r * n_proc_cols + c

            col_end = col_cumsum[c]
            matrix[row_start : row_end , col_start : col_end] = A.reshape((imax_local_array[index], jmax_local_array[index]))
            
            plt.plot([col_start, col_end], [row_end, row_end], 'b-', linewidth = 1)
            plt.plot([col_end, col_end], [row_start, row_end], 'b-', linewidth = 1)
            
            col_start = col_end
        row_start = row_end   

    matrix_store[it, :, :] = matrix

plt.title("Domain Decomposition for %i x %i environment" % (imax_global, jmax_global))
plt.show()

def create_animation():
    fig = plt.figure()
    plt.title('Game of Life')
    plt.xlabel('i')
    plt.ylabel('j')
    imgs = []
    
    for it in range(iter):
        img = plt.imshow(matrix_store[it], 'binary', animated=True)
        title = plt.title("%i generations" % it, animated=True)
        imgs.append([img, title])    
    ani = anim.ArtistAnimation(fig, imgs, interval=50, blit=True)
    ani.save('game.mp4')
    plt.close(fig)  # prevent final frame plot from showing up inline below

    return ani

print("Computational time for Game of Life = %f" % np.mean(time_array))
ani = create_animation()
print('Preparing mp4 gif (please wait a moment.)')
ani.save('game.mp4')