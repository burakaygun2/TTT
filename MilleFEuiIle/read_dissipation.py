import h5py

with h5py.File("../results/output.h5", "r") as f:
    dissipation = f["Periodic/Dissipation/average_dissipation_layer_5"][()]
    print(dissipation)