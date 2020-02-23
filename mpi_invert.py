import hazel
import h5py
import sys
conf_file = sys.argv[1]

iterator = hazel.Iterator(use_mpi=True)
mod = hazel.Model(conf_file, rank=iterator.get_rank(),working_mode='inversion')
iterator.use_model(model=mod)
iterator.run_all_pixels()
