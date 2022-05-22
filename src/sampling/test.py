import numpy as np
import matplotlib.pyplot as plt
import copy
# From https://github.com/mrphysics-bonn/skipped-caipi:
from skippedcaipi.skippedcaipi import elementary_sampling, get_trajectory_indices, plot_parabola_connection, get_zblips, get_zblipcycle

import skippedcaipi_sampling

Ry=1 # Undersampling factor along y (primary phase encode direction, w.l.o.g.)
Rz=6 # Undersampling factor along z (slice direction, w.l.o.g.)
Dz=2 # CAIPI shift along z
S=1  # Segmentation factor (typically 1, for fMRI)

matrix_size_y = 64
matrix_size_z = 30

# Create the blipped-CAIPI instance linked to the CAIPI pattern (Segmentation Factor = 1)
blippedcaipi = skippedcaipi_sampling.skippedcaipi_sampling(matrix_size_y, matrix_size_z, Ry, Rz, Dz, SegmentationFactor=S)

plt.figure(figsize=(5,10))
blippedcaipi.plot()
plt.title('The blipped-CAIPI instance')
plt.show()
