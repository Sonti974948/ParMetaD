import os
import numpy as np
import h5py

# Set your working directory
path = '/expanse/lustre/scratch/ssonti/temp_project/amber_learn/chignolin_tutorial/westpa_tutorials/tutorial7.3-chignolin/ParGaMD/ParMetaD/run_2'
os.chdir(path)

h5_path = "west.h5"

# Initialize lists to hold the columns
all_cv = []
all_rg = []
all_rbias = []
all_weights = []

print("Extracting data from west.h5 and COLVAR files...")

# Open the HDF5 file
with h5py.File(h5_path, "r") as h5:
    
    # Get total number of iterations (-1 because index 0 is initial state)
    num_iterations = len(h5["summary"]) - 1
    
    for iter_idx in range(1, num_iterations + 1):
        
        # Access the specific iteration in west.h5 to get the weights
        iter_group = h5[f"iterations/iter_{iter_idx:08d}"]
        weights = iter_group["seg_index"]["weight"]
        num_walkers = len(weights)
        
        for seg_idx in range(num_walkers):
            
            # 1. Get the weight for this specific walker
            walker_weight = weights[seg_idx]
            
            # 2. Construct the path to the COLVAR file
            # Format: traj_segs/000001/000000/COLVAR
            colvar_path = f"traj_segs/{iter_idx:06d}/{seg_idx:06d}/COLVAR"
            
            try:
                # Load COLVAR file. np.loadtxt automatically ignores the '#! FIELDS' header
                data = np.loadtxt(colvar_path)
                
                # Handle cases where COLVAR might only have 1 frame
                if data.ndim == 1:
                    data = data.reshape(1, -1)
                
                # In your previous message: 
                # Col 0=time, Col 1=cv, Col 2=bias, Col 3=rct, Col 4=rbias
                cv = data[:, 1]
                rg = data[:,2]
                rbias = data[:, 3]
                
                
                # Create an array that repeats the walker weight for every frame (e.g., 500 times)
                num_frames = len(cv)
                weight_array = np.full(num_frames, walker_weight)
                
                # Append to master lists
                all_cv.extend(cv)
                all_rg.extend(rg)
                all_rbias.extend(rbias)
                all_weights.extend(weight_array)
                
            except Exception as e:
                # If a COLVAR file is missing or crashed, gracefully skip it
                continue
                
        # Optional: Print progress every 10 iterations to ensure it's not frozen
        if iter_idx % 10 == 0:
            print(f"Processed Iteration {iter_idx}/{num_iterations}")

# --- Formatting and Saving ---
print("Aggregating master array...")

# Column stack groups the lists into a 2D array (Row = Frame, Cols = CV, rbias, weight)
master_data = np.column_stack((all_cv,all_rg, all_rbias, all_weights))

print(f"Total frames extracted: {len(master_data)}")
print("Saving to weights_master.dat...")

# Save to text file. Using scientific notation (%e) to preserve the tiny walker weights
np.savetxt(
    'weights_master.dat', 
    master_data, 
    header="CV RG metad.rbias walker_weight", 
    fmt='%.6e %.6e %.6e %.10e'
)

print("Done! weights_master.dat has been created.")

