#!/usr/bin/env python
#-*- coding:utf-8 -*-

"""skippedcaipi_sampling.py: This class provides basic functionality to generate a CAIPIRINHA sampling mask and a skipped-CAIPI [1] view ordering list from the specified phase encoding matrix size and CAIPI parameters. This can be used for blipped-CAIPI [2], shot-selective CAIPI [3,4] (without z-blips) or any other skipped-CAIPI [1] sampling.

Note: skipped-CAIPI was developed for 3D-EPI, but in principle it should be applicable to SMS-EP as well [5]. The blips should remain the same, only the initial phase encoders must be realized differently in SMS-EPI (via RF phase?).

1. Stirnberg, R., & Stöcker, T. (2021). Segmented K‐space blipped‐controlled aliasing in parallel imaging for high spatiotemporal resolution EPI. Magnetic Resonance in Medicine, 85(3), 1540–1551. https://doi.org/10.1002/mrm.28486

2. Setsompop, K., Gagoski, B. A., Polimeni, J. R., Witzel, T., Wedeen, V. J., & Wald, L. L. (2012). Blipped-controlled aliasing in parallel imaging for simultaneous multislice echo planar imaging with reduced g-factor penalty. Magnetic Resonance in Medicine, 67(5), 1210–1224. https://doi.org/10.1002/mrm.23097

3. Poser, B. A., Ivanov, D., Kannengiesser, S. A., Uludağ, K., & Barth, M. (2014). Accelerated 3D EPI using 2D blipped-CAIPI for high temporal and/or spatial resolution. Proceedings of the International Society of Magnetic Resonance in Medicine, 22, 1506.

4. Hendriks, A. D., D’Agata, F., Raimondo, L., Schakel, T., Geerts, L., Luijten, P. R., Klomp, D. W. J., & Petridou, N. (2020). Pushing functional MRI spatial and temporal resolution further: High-density receive arrays combined with shot-selective 2D CAIPIRINHA for 3D echo-planar imaging at 7 T. NMR in Biomedicine, 33(5), 1–13. https://doi.org/10.1002/nbm.4281

5. Zahneisen, B., Poser, B. A., Ernst, T., & Stenger, V. A. (2014). Three-dimensional Fourier encoding of simultaneously excited slices: Generalized acquisition and reconstruction framework. Magnetic Resonance in Medicine, 71(6), 2071–2081. https://doi.org/10.1002/mrm.24875

"""

__author__ = "Rüdiger Stirnberg"
__maintainer__ = "Rüdiger Stirnberg"
__email__ = "ruediger.stirnberg@dzne.de"
__version__ = "1.0"

import numpy as np
import matplotlib.pyplot as plt
import copy

# From https://github.com/mrphysics-bonn/skipped-caipi:
from skippedcaipi import elementary_sampling, get_trajectory_indices, plot_parabola_connection, get_zblips, get_zblipcycle

# For calling from linux shell
# Example usage:  
#   $ python skippedcaipi_sampling.py 1 6 2 1 64 30
if __name__ == "__main__":
    import sys
    import skippedcaipi_sampling

    Ry = int(sys.argv[1]) # Undersampling factor along y (primary phase encode direction, w.l.o.g.)
    Rz = int(sys.argv[2]) # Undersampling factor along z (slice direction, w.l.o.g.)
    Dz = int(sys.argv[3]) # CAIPI shift along z
    S = int(sys.argv[4])  # Segmentation factor (typically 1, for fMRI)

    matrix_size_y = int(sys.argv[5])
    matrix_size_z = int(sys.argv[6])

    # Create an instance of the blipped-CAIPI object 
    blippedcaipi = skippedcaipi_sampling.skippedcaipi_sampling(matrix_size_y, matrix_size_z, Ry, Rz, Dz, SegmentationFactor=S)

    #for echo in range(blippedcaipi.epi_factor(0)//2):
    #    print(echo, blippedcaipi.indices[0][echo,:])

    # write ky/kz sampling pattern to .mat file
    blippedcaipi.tomatfile()

class skippedcaipi_sampling:
    def __init__(self, matrix_size_y, matrix_size_z, Ry, Rz, CaipiShiftZ, SegmentationFactor=1):
        self.matrix_size = [matrix_size_y, matrix_size_z]
        self.R = [Ry,Rz]
        self.D = CaipiShiftZ
        self.S = SegmentationFactor
        
        self.mask = self.samplingmask()
        self.samples = np.sum(self.mask)
        
        # update actual matrix size
        self.matrix_size = self.mask.shape
        
        self.sampling_repeats_z = self.matrix_size[1]//self.R[1]
        
        # Update sampling according to selected segmentation factor
        self.update_sampling()
        
    def update_sampling(self, SegmentationFactor=None):
        if SegmentationFactor is not None:
            self.S = SegmentationFactor
        
        self.shots_to_measure = self.sampling_repeats_z * self.S
            
        self.indices = self.viewordering()
        
        self.zblips = get_zblips(self.R[0], self.R[1], self.D, self.S)
        self.zblip_cycle = get_zblipcycle(self.R[1],np.min(self.zblips))
        
    def update_shotselective(self):
        self.update_sampling(SegmentationFactor=1)
        blippedcaipi_zblip_cycle = self.zblip_cycle
        
        self.update_sampling(SegmentationFactor=blippedcaipi_zblip_cycle)

    def has_zblips(self):
        return (self.zblip_cycle!=1)
    
    def epi_factor(self, shot=0):
        if shot < 0 or shot >= self.shots_to_measure:
            return -1 # in valid
        else:
            return self.indices[shot].shape[0]
        
    def samplingmask(self):
        Rtot = np.prod(self.R)
        sampling_cell_repeats_y = int(np.ceil(self.matrix_size[0]*1.0/Rtot))
        sampling_cell_repeats_z = int(np.ceil(self.matrix_size[1]*1.0/Rtot)) 

        # Along y, it is not necessary that the matrix size is divisable by Rtot.
        mask = elementary_sampling(self.R[0], self.R[1], self.D, sampling_cell_repeats_y-1).T
        mask = mask[:self.matrix_size[0],:]
        
        # Along z, it can make sense that the matrix size (number of "slices") are divisable by Rz.
        mask = np.tile(mask, (1,sampling_cell_repeats_z))
        mask = mask[:,:int(np.ceil(self.matrix_size[1]*1.0/self.R[1])*self.R[1])]
    
        return mask
    
    def viewordering(self):            
        
        # First S shots ...
        sampling_template = [get_trajectory_indices(self.mask.T, segmentation=self.S, shot=s) for s in range(self.S)]
        sampling_indices = []
        # ... that are repeated for subsequent k-space partitions
        for instance in range(self.sampling_repeats_z):
            sampling_indices += copy.deepcopy(sampling_template)
            partition_offset = self.R[1] * instance
            for s in range(self.S):
                # Add partition offset
                sampling_indices[-self.S+s][:,0] += partition_offset
                
        return sampling_indices
    
    def plot(self, axes=None, shots=None):
        if not axes:
            axes = plt.gca()
            
        if shots==None:
            shots = range(self.shots_to_measure)

        axes.imshow(self.mask.T, cmap='gray')

        # color of trajectories encodes time:
        colors = plt.cm.hot(np.linspace(0,1,self.samples-1))
        sample = 0
        
        for s in shots:
            if s>=0 and s < self.shots_to_measure:
                actual_epi_factor = self.indices[s].shape[0]
                for echo in range(actual_epi_factor-1):
                    if not self.has_zblips():
                        plt.plot(self.indices[s][echo:echo+2,1], self.indices[s][echo:echo+2,0], color=colors[sample])
                    else:
                        #A little more realistic phase encode trajectory between echoes is plotted:
                        plot_parabola_connection(self.indices[s][echo,::-1], self.indices[s][echo+1,::-1], axis=axes, num=11, color=colors[sample], bRotate=False)
                    sample += 1

        axes.set_ylabel('Parts')
        axes.set_xlabel('Lines')

    def tomatfile(self):
        import scipy.io as spio
        mdict = {"indices": self.indices[0], "matrix_size": self.matrix_size, "R": self.R,
            "D": self.D, "S": self.S, "mask": self.mask, "samples": self.samples, 
            "sampling_repeats_z": self.sampling_repeats_z}
        spio.savemat("caipi.mat", mdict)

