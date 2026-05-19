import h5py

def count_walkers(west_path="west.h5"):
    total_walkers = 0

    with h5py.File(west_path, "r") as f:
        if "iterations" not in f:
            print("No 'iterations' group found in west.h5")
            return

        iters_group = f["iterations"]

        print("Walkers per iteration:\n")

        for iter_key in sorted(iters_group.keys()):
            iter_group = iters_group[iter_key]

            if "seg_index" not in iter_group:
                continue

            seg_index = iter_group["seg_index"]
            n_walkers = len(seg_index)

            print(f"{iter_key}: {n_walkers}")
            total_walkers += n_walkers

    print("\nTotal walkers across all iterations:", total_walkers)


if __name__ == "__main__":
    count_walkers("west.h5")

