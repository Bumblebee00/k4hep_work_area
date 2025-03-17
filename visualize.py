import sys
import uproot
import awkward as ak
import matplotlib.pyplot as plt

if len(sys.argv) != 2:
    print("Usage: python3 visualize.py <edm4hep root file>")
    exit()

# Open the file and get the TTree of events
file = uproot.open(sys.argv[1])
tree = file["events"]

# Extract the hit energies from the ECalBarrelCollection branch
energy = tree["ECalBarrelCollection/ECalBarrelCollection.energy"].array()
flat_energy = ak.flatten(energy)  # create a 1D array

mean = ak.mean(flat_energy)
std = ak.std(flat_energy)
print(f"energy mean: {mean*1e3:.6f} MeV, standard deviation: {std:.6f}")

print("Number of hits caused by the first, second, etc... event:")
print(ak.num(energy))  # Show the number of hits in each event
print("Number of events: ", len(energy))
print("Total number of hits: ", len(flat_energy))

# Figure 1: Energy histogram
fig1 = plt.figure(figsize=(10, 6))
plt.hist(flat_energy, bins=100, histtype="step", color="blue")
plt.xlabel("Hit Energy [GeV]")
plt.ylabel("Number of Hits")
plt.title("ECalBarrel Hit Energy Distribution")
plt.show()
