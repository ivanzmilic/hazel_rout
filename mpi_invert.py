import hazel
import h5py

iterator = hazel.Iterator(use_mpi=True)
mod = hazel.Model('conf_multi.ini', rank=iterator.get_rank(),working_mode='inversion')
iterator.use_model(model=mod)
iterator.run_all_pixels()
