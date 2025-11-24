"""Author: Andrew Martin
Creation date: 04/10/2025

Script to handle computation of Mutual Information for the Vertical Cloud Fraction datasets using the following methods:
    + Holmes estimator of full (2x50)-dimensional dataset

Then, a null-hypothesis distribution is generated from independent samples, and a p value computed that the computed MI value is significant (not independent data).
Then, the std of the MI estimator is estimated using the method outlined in Holmes et al. (2019).
"""

from ..definitions import indices
from ..holmes_et_al_2019 import python_interface as holmes

from scipy import stats 

from dataclasses import dataclass
import xarray as xr
import numpy as np
import os



SITES = ("ny-alesund", "hyytiala", "juelich", "munich",)



@dataclass(kw_only=True, frozen=True)
class Args:
    dir_vcfs: str
    site: str
    index_function_name: str
    array_index: int | None
    fpath_out: str
    K: int
    n_bootstraps: int
    n_splits: int
    n_B_repeats: int
    seed: int


def parse_args() -> Args:
    import argparse
    parser = argparse.ArgumentParser(description="Collocate ATL09 and CloudNet data")

    # Required arguments (no default values in dataclass)
    parser.add_argument(
        "--dir-vcfs",
        type=str,
        required=True,
        help="Directory path for the vcfs-per-event netcdf files"
    )

    parser.add_argument(
        "--site",
        type=str,
        choices = SITES,
        required=True,
        help="Site identifier"
    )

    parser.add_argument(
        "--dir-out",
        type=str,
        required=True,
        help="Directory path to location for saving output Mutual-Information-cube netcdf files"
    )

    parser.add_argument(
        "--index-function",
        choices = indices.INDEX_FUNCTIONS.keys(),
        required=True,
        help="Name of the function mapping job array indices to (R,tau) pairs"
    )

    parser.add_argument(
        '--job-array-index',
        type=int,
        required=False,
        help='Job array index -- only runs a single computation'
    )
    parser.add_argument(
        "-K",
        type=int,
        required=True,
        help="The order of the Holmes estimator (finds the kth nearest neighbour)"
    )

    parser.add_argument(
        "--n-bootstraps",
        type=int,
        required=True,
        help="The number of bootstrapped independent samples to test the significance of the computed MI values."
    )
    parser.add_argument(
        "--n-splits",
        type=int,
        required=True,
        help="The number of non-overlapping splits made to the VCF pairs to estimate a standard deviation of the MI value."
    )
    parser.add_argument(
        "--n-B-repeats",
        type=int,
        required=True,
        help="The number of times variance of the MI estimator should be computed, to obtain the maximum likelihood estimate of B."
    )

    parser.add_argument(
        "--seed",
        type=int,
        required=False,
        help="Seed to initialise random number generator for bootstrapping"
    )

    parsed_args = parser.parse_args()
    
    assert os.path.isdir( dir_vcfs := parsed_args.dir_vcfs)
    assert os.path.isdir( dir_out := parsed_args.dir_out )

    site = parsed_args.site
    K = parsed_args.K

    fname = f"MI_{site}_K{K}"
    if (job_array_index:=parsed_args.job_array_index) is not None:
        fname += f"_job{job_array_index}"  
    fname += ".nc"
    fpath_out = os.path.join(
        dir_out,
        fname
    )

    return Args(
        dir_vcfs = dir_vcfs,
        site = site,
        fpath_out = fpath_out,
        index_function_name = parsed_args.index_function,
        array_index = job_array_index, 
        K = K,
        n_bootstraps = parsed_args.n_bootstraps,
        n_splits = parsed_args.n_splits,
        n_B_repeats = parsed_args.n_B_repeats,
        seed = parsed_args.seed,
    )



def iterate_index_func(index_func):
    """Generator that stops when InvalidIndexError is raised"""
    i = 0
    while True:
        try:
            yield index_func(i)
            i += 1
        except indices.InvalidIndexError:
            return  # Stops the generator



def parameterisation_to_netcdf_fname(param: indices.Parametrisation, site: str) -> str:
    """takes a Parametrisation output of an index function and converts it to a filename according to the outputs from compute_vcfs_per_event.py
    """
    R_km = param.distance_km
    tau_s = int( param.tau.total_seconds() )
    return f"vcfs-per-event_{site}_{R_km:3.3f}km_{tau_s:06}s.nc"



def generate_independent_MI_values(X: np.ndarray, Y: np.ndarray, n_bootstraps: int, K: int, rng: np.random.Generator) -> np.ndarray:
    """Function that given data arrays of shape (d,N_events), bootstraps independent data of shape (d,N_events) nb_bootstraps times, and computes the MI each time, generating a null-hypothesis-distribution that the original data are indepenedent VCF profiles.
    """
    assert X.shape[1] == Y.shape[1]
    MI_independent_values = list()
    N_events = X.shape[1]
    for _ in range(n_bootstraps):
        # generate two sets of integers in [0, N_events) that are strictly elementwise not equal
        randintX = np.zeros(N_events, dtype=int)
        randintY = np.zeros(N_events, dtype=int)
        while (randintX == randintY).any():
            update_mask = (randintX == randintY)
            n_update = update_mask.sum()
            randintX[update_mask] = rng.integers(low=0, high=N_events, size=n_update)
            randintY[update_mask] = rng.integers(low=0, high=N_events, size=n_update)
        
        X_bootstrap = X[:, randintX]
        Y_bootstrap = Y[:, randintY]
        
        MI_independent_values.append(
            holmes.call_MI_xnyn(
                X = X_bootstrap,
                Y = Y_bootstrap,
                K = K
            )
        )
    return np.array(MI_independent_values)



def estimate_B_ML(X: np.ndarray, Y: np.ndarray, n_splits: int, n_repeats: int, K: int, rng: np.random.Generator) -> [float, float]:
    """Function that computes the Maximum Likelihood B value for modelling the MI estimate variance, using both the literature equation and the correction I derived.
    This works by creating n_splits non-overlapping partitions of the data and computing a variance for the MI estimates.
    This is performed n_repeats times, and a ML value of B is computed.
    """
    assert X.shape[1] == Y.shape[1]
    N_events = X.shape[1]

    variance_estimates = list()
    n_i_values = list()
    for _ in range(n_repeats):
        shuffled_indices = rng.permutation(np.arange(N_events))
        non_overlapping_indices = np.array_split(shuffled_indices, n_splits)

        MI_subsets = list()
        for indices in non_overlapping_indices:
            X_subset = X[:, indices]
            Y_subset = Y[:, indices]
            MI_subsets.append(
                holmes.call_MI_xnyn(
                    X = X_subset,
                    Y = Y_subset,
                    K = K
                )
            )
        variance_estimates.append(np.var(MI_subsets))
        n_i_values.append(n_splits)

    B_numerator = N_events * np.sum([
        (n_i - 1) / n_i * var
        for var, n_i in zip(variance_estimates, n_i_values)
    ])
    B_denominator_lit = np.sum([ n_i - 1 for n_i in n_i_values ])
    B_denominator_corrected = np.sum([ (n_i - 3)/2 for n_i in n_i_values ])
    return float(B_numerator / B_denominator_lit), float(B_numerator / B_denominator_corrected)



def process_single_parametrisation(
    params: indices.Parametrisation,
    args: Args,
    RNG: np.random.Generator,
    verbose: bool = True,
) -> xr.Dataset:
    fpath = os.path.join(
        args.dir_vcfs,
        parameterisation_to_netcdf_fname(params, args.site)
    )
    if verbose: 
        print("")
        print(f"LOADING {fpath=}", end=", ")
    ds = xr.load_dataset(fpath).transpose("height", "collocation_event")
    if verbose: print(f"loading-success")

    R_km = params.distance_km
    tau_s = int( params.tau.total_seconds() )
    
    N_events = int(ds.collocation_event.size)
    # data required for Holmes estimators
    X = ds.vcf_atl09.data
    Y = ds.vcf_cloudnet.data

    data_vars = dict()

    # use my holmes implementation of MIEstimate
    MI_estimate = holmes.MIEstimate.from_XYKMn_with_RNG(
        X = X,
        Y = Y,
        K = args.K,
        M = args.n_B_repeats, # number of repeated computations of sigma_KSG_i for computing sigma_KSG(n_samples)
        n_splits = args.n_splits, # number of splits, n_i, per computation of sigma_KSG_i
        n_samples = N_events, # assertion to ensure data passed to estimators correctly
        RNG = RNG, # controls M permutations of X and Y for computing std_MI
    )

    # compute a distribution of MI values for independent shufflings of the input data
    MI_independent = generate_independent_MI_values(
        X = X,
        Y = Y,
        n_bootstraps = args.n_bootstraps,
        K = args.K,
        rng = RNG
    )

    # compute p-value using an unpaired Welch's t-test (non-equal variance) to show that the computed MI values are not drawn from the same distribution as the independent MI samples.
    pvalue_independent = stats.ttest_ind_from_stats(
        mean1 = MI_estimate.MI,
        std1 = MI_estimate.std,
        nobs1 = MI_estimate.ddof,
        mean2 = np.mean(MI_independent),
        std2 = np.std(MI_independent),
        nobs2 = args.n_bootstraps,
        equal_var=False,
    ).pvalue

    if verbose: 
        print(f"I_KSG({R_km} km, {tau_s} s) = {MIEstimate.MI}")
        print(f"std_KSG =  {MI_estimate.std}")
        print(f"independent p-value=  {pvalue_independent}")
    
    data_vars["MI"] = MI_estimate.MI
    data_vars["std"] = MI_estimate.std
    data_vars["pvalue_independent"] = pvalue_independent
    # including n_events and n_profiles
    data_vars["N_events"] = N_events
    data_vars["N_profiles"] = int(
        (ds.n_profiles_atl09 * ds.n_profiles_cloudnet).sum()
    )
    data_vars["n_splits_std"] = MI_estimate.n_splits
    data_vars["std_ddof"] = MI_estimate.ddof

    # create results dataset
    if verbose: print(f"creating results dataset: ", end="")
    results_dataset = xr.Dataset(
        data_vars = data_vars,
        coords = {
            "height": ds.coords["height"]
        }
    ).expand_dims(
        dim = {
            "R_km": [R_km],
            "tau_s": [tau_s],
            "K": [args.K]
        }
    )
    if verbose: print("success")
    return results_dataset




def main(args: Args):
    """For a given site, loads all of the vcfs per event datasets, and calculates the 4 MI values for each.
    Then, combines all datasets (including N_events and N_profiles) into a dataset with dimensions (R_km, tau_s, height).
    """
    fpath_out = args.fpath_out 
    print(f"{fpath_out=}")

    random_seed = args.seed if args.seed is not None else np.random.default_rng().integers(low=0, high=np.iinfo(np.int64).max)
    #random_seed = np.random.default_rng().integers(low=0, high=np.iinfo(np.int64).max)
    print(f"RANDOM SEED: {random_seed}")
    RNG = np.random.default_rng(random_seed)

    index_function = indices.INDEX_FUNCTIONS[args.index_function_name]

    if (job_index := args.array_index) is not None:
        MI_ds = process_single_parametrisation(
            params = index_function(job_index),
            args = args,
            RNG=RNG,
        )
    else:
        MI_datasets = list()
        for n, params in enumerate(iterate_index_func(index_function)):
            results_dataset = process_single_parametrisation(
                params = params,
                args = args,
                rng = RNG,
            )
            MI_datasets.append(results_dataset)

        print("\n"*3)
        # combine all individual results datasets by their R and tau coordinates
        print("COMBINING datasets: ", end="")
        MI_ds = xr.combine_by_coords(MI_datasets)
        print("success")

    print(f"SAVING dataset to {fpath_out=}")
    MI_ds.to_netcdf(fpath_out)
    print("SAVING SUCCESS")



if __name__ == "__main__":
    args = parse_args()
    print(args)
    main(args)
    print("SCRIPT SUCCESS")
