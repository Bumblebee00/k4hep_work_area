import sys
import uproot
import awkward as ak
import matplotlib.pyplot as plt
import numpy as np
from sklearn.cluster import KMeans, DBSCAN, MeanShift


if len(sys.argv) != 2:
    print("Usage: python3 cluster.py <edm4hep root file>")
    exit()

# Open the file and get the TTree of events
file = uproot.open(sys.argv[1])
tree = file["events"]

# Extract the hit energies from the ECalBarrelCollection branch
energy = tree["ECalBarrelCollection/ECalBarrelCollection.energy"].array()
positions_x = tree["ECalBarrelCollection/ECalBarrelCollection.position.x"].array()
positions_y = tree["ECalBarrelCollection/ECalBarrelCollection.position.y"].array()
positions_z = tree["ECalBarrelCollection/ECalBarrelCollection.position.z"].array()

flat_energy = ak.flatten(energy)  # create a 1D array
flat_x = ak.flatten(positions_x)
flat_y = ak.flatten(positions_y)
flat_z = ak.flatten(positions_z)

print("Number of hits caused by the first, second, etc... event:")
print(ak.num(energy))  # Show the number of hits in each event
 
# does nothing, just to plot the original
def identity(x,y,z,e):
    newposs = []
    newes = []

    for i in range(len(x)):
        newposs.append((x[i], y[i], z[i]))
        newes.append(e[i])
    return newposs, newes, "Original"

# divide the space in cubes of size 10, and for each cube, calculate the baricenter of the hits and the total energy
def grid(x,y,z,e):
    newposs = []
    newes = []

    grid = {}

    for i in range(len(x)):
        key = (int(x[i]/10), int(y[i]/10), int(z[i]/10))
        if key not in grid:
            grid[key] = ([], [])
        grid[key][0].append((x[i], y[i], z[i]))
        grid[key][1].append(e[i])
    
    for key in grid:
        newposs.append((np.mean([x[0] for x in grid[key][0]]), np.mean([x[1] for x in grid[key][0]]), np.mean([x[2] for x in grid[key][0]])))
        newes.append(np.sum(grid[key][1]))

    return newposs, newes, "Grouped by 10x10x10 cubes"

def grid20(x,y,z,e):
    newposs = []
    newes = []

    grid = {}

    for i in range(len(x)):
        key = (int(x[i]/20), int(y[i]/20), int(z[i]/20))
        if key not in grid:
            grid[key] = ([], [])
        grid[key][0].append((x[i], y[i], z[i]))
        grid[key][1].append(e[i])
    
    for key in grid:
        newposs.append((np.mean([x[0] for x in grid[key][0]]), np.mean([x[1] for x in grid[key][0]]), np.mean([x[2] for x in grid[key][0]])))
        newes.append(np.sum(grid[key][1]))

    return newposs, newes, "Grouped by 20x20x20 cubes"


# wrapper for the sklearn DBSCAN algorithm
# fail
def mydbscan(x, y, z, e):
    # shuffle x,y,z by +- 0.1
    # x = flat_x + np.random.uniform(-0.1, 0.1, len(flat_x))
    # y = flat_y + np.random.uniform(-0.1, 0.1, len(flat_y))
    # z = flat_z + np.random.uniform(-0.1, 0.1, len(flat_z))
    # e = flat_energy

    newposs = []
    newes = []

    model = DBSCAN(eps=30, min_samples=1)
    points = np.column_stack((x, y, z))
    model.fit(points)
    labels = model.labels_
    n_clusters = len(set(labels)) - (1 if -1 in labels else 0)
    # calculate centroid of each cluster
    for i in range(n_clusters):
        mask = (labels == i)
        if np.any(mask):
            newposs.append(np.mean(points[mask], axis=0))
            newes.append(np.sum(np.array(e)[mask]))
    return newposs, newes, f"DBSCAN clustering (eps=30, min_samples=1)"

def kmeans4d(x, y, z, e):
    points = np.column_stack((x, y, z, e))
    k = KMeans(n_clusters=150, n_init='auto', init='random').fit(points)
    newposs = k.cluster_centers_[:, :3]
    newes = k.cluster_centers_[:, 3]
    return newposs, newes, "K-means clustering 4D"

def kmeans(x, y, z, e):
    points = np.column_stack((x, y, z))
    k =KMeans(n_clusters=150, n_init='auto', init='random').fit(points)
    newposs = k.cluster_centers_
    newes = np.zeros(150)
    for i in range(len(points)):   
        newes[k.labels_[i]] += e[i]
    return newposs, newes, "K-means clustering"

def meanshift(x, y, z, e):
    points = np.column_stack((x, y, z))
    ms = MeanShift(bandwidth=30).fit(points)
    newposs = ms.cluster_centers_
    newes = np.zeros(len(newposs))
    for i in range(len(points)):
        newes[ms.labels_[i]] += e[i]
    return newposs, newes, "MeanShift clustering (sklearn)"

# takes a array of position and energy of hits, and creates a pointcloud of smaller dimension that preserves key calorimetric observables
strategies = [identity, grid, kmeans, kmeans4d, meanshift, mydbscan, grid20]

print(f"{'Strategy':<30} | {'Points':>6} | {'Energy Mean (GeV)':>17} | {'Std Dev':>10} | {'RMS':>10} | {'Baricenter X':>15} | {'Baricenter Y':>15} | {'Baricenter Z':>15} | {'Total Energy':>20}")
print("-" * 160)


for f in strategies:
    newposs, newes, name = f(np.array(positions_x[0]), np.array(positions_y[0]), np.array(positions_z[0]), energy[0])
    energy_mean = np.mean(newes)
    energy_std = np.std(newes)
    # weighted baricenter by the energy of each hit
    baricenter = (np.sum([x[0]*e for x, e in zip(newposs, newes)])/np.sum(newes), np.sum([x[1]*e for x, e in zip(newposs, newes)])/np.sum(newes), np.sum([x[2]*e for x, e in zip(newposs, newes)])/np.sum(newes))
    # caluclate root mean square radius
    rms = np.sqrt(np.sum([np.linalg.norm(np.array(x)-np.array(baricenter))**2*e for x, e in zip(newposs, newes)])/np.sum(newes))
    
    print(f"{name:<30} | {len(newposs):>6} | {energy_mean:>17.5f} | {energy_std:>10.5f} | {rms:>10.5f} | "
          f"{baricenter[0]:>15.6f} | {baricenter[1]:>15.6f} | {baricenter[2]:>15.6f} | {sum(newes):>20.10f}")


    # create a 3d matplotlib scatter plot of the new data, with equal axis scale
    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection='3d')
    scatter = ax.scatter([x[0] for x in newposs], [x[1] for x in newposs], [x[2] for x in newposs], c=newes, s=2, alpha=0.6, 
                        cmap='brg', vmin=0, vmax=0.0025)
    plt.colorbar(scatter, ax=ax, label='Energy [GeV]')
    ax.set_xlabel('X Position')
    ax.set_ylabel('Y Position')
    ax.set_zlabel('Z Position')
    ax.set_title(f'3D Distribution of Detector Hits for one run ({name})')
    x_limits = ax.get_xlim()
    y_limits = ax.get_ylim()
    z_limits = ax.get_zlim()
    x_range = abs(x_limits[1] - x_limits[0])
    y_range = abs(y_limits[1] - y_limits[0])
    z_range = abs(z_limits[1] - z_limits[0])
    max_range = max(x_range, y_range, z_range)
    x_center = (x_limits[0] + x_limits[1]) / 2
    y_center = (y_limits[0] + y_limits[1]) / 2
    z_center = (z_limits[0] + z_limits[1]) / 2
    ax.set_xlim(x_center - max_range/2, x_center + max_range/2)
    ax.set_ylim(y_center - max_range/2, y_center + max_range/2)
    ax.set_zlim(z_center - max_range/2, z_center + max_range/2)


plt.show()