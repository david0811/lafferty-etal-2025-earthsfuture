import json

import jax.numpy as jnp
import numpy as np
import pandas as pd
from SALib import ProblemSpec
from scipy.stats import gaussian_kde, rankdata
from utils.global_paths import project_data_path

from water_balance_jax import (
    construct_Kpet_crop,
    construct_Kpet_gen,
    wbm_jax,
)

def calculate_delta_full(Y, X):
    """Calculate delta sensitivity indices using the full dataset without random sampling.
    For original implementation, see: https://salib.readthedocs.io/en/latest/_modules/SALib/analyze/delta.html
    
    Parameters:
    -----------
    Y : numpy.ndarray
        Model outputs, shape (N,)
    X : numpy.ndarray
        Model inputs, shape (N, D) where D is number of parameters
        
    Returns:
    --------
    delta : numpy.ndarray
        Delta sensitivity indices for each parameter
    """
    N = len(Y)
    D = X.shape[1]
    delta = np.zeros(D)
    
    # Create grid for PDF estimation
    Ygrid = np.linspace(np.min(Y), np.max(Y), 100)
    
    # Estimate the unconditional PDF
    fy = gaussian_kde(Y, bw_method="silverman")(Ygrid)
    
    # Equal frequency partition calculation (similar to original code)
    exp = 2.0 / (7.0 + np.tanh((1500.0 - N) / 500.0))
    M = int(np.round(min(int(np.ceil(N**exp)), 48)))
    m = np.linspace(0, N, M + 1)
    
    for i in range(D):
        X_i = X[:, i]
        xr = rankdata(X_i, method="ordinal")
        
        d_hat = 0.0
        for j in range(len(m) - 1):
            ix = np.where((xr > m[j]) & (xr <= m[j + 1]))[0]
            nm = len(ix)
            
            Y_ix = Y[ix]
            if np.ptp(Y_ix) != 0.0:
                fyc = gaussian_kde(Y_ix, bw_method="silverman")(Ygrid)
                fy_ = np.abs(fy - fyc)
            else:
                fy_ = np.abs(fy)
            
            d_hat += (nm / (2 * N)) * np.trapezoid(fy_, Ygrid)
        
        delta[i] = d_hat
    
    return delta

def wbm_sobol(
    ix,
    iy,
    forcing,
    eval,
    tas_delta,
    prcp_factor,
    Kpet_name,
    experiment_name,
    N,
    store_sims,
    save_name,
):
    """
    Perform a single gridpoint Sobol SA using parameter file in `experiment_name`.
        - ix, iy assuming eCONUS (whole domain)
        - Metrics include mean, SD, range, and RMSE against obs if `eval` = NLDAS, SMAP
        - CC can be imposed by `tas_delta` and `prcp_factor`
    """
    ##################
    # Get forcing
    ##################
    # Read data
    if forcing == "SMAP":
        forcing_data = np.load(
            f"{project_data_path}/WBM/calibration/eCONUS/{forcing}/inputs.npz"
        )
    else:
        forcing_data = np.load(
            f"{project_data_path}/WBM/calibration/eCONUS/VIC/inputs.npz"
        )  # identical for SA purposes

    # Select gridpoint
    tas_in = forcing_data["tas"][ix, iy, :]
    prcp_in = forcing_data["prcp"][ix, iy, :]
    lai_in = forcing_data["lai"][ix, iy, :]
    phi = forcing_data["lats"][iy]

    ###################
    # Get observations
    ###################
    if eval == "SMAP":
        obs = np.load(
            f"{project_data_path}/WBM/calibration/eCONUS/{forcing}/{forcing}_validation.npy"
        )
        obs = obs[ix, iy, :]
        obs_centered = obs - jnp.mean(obs)
    elif eval == "NLDAS":
        obs_list = ["VIC", "NOAH", "MOSAIC"]
        obs = [
            np.load(
                f"{project_data_path}/WBM/calibration/eCONUS/{obs}/{obs}_validation.npy"
            )
            for obs in obs_list
        ]
        obs = [obs_tmp[ix, iy, :] for obs_tmp in obs]
        obs_centered = [obs_tmp - jnp.mean(obs_tmp) for obs_tmp in obs]

    ################
    # Get params
    ################
    params = np.loadtxt(
        f"{project_data_path}/WBM/SA/{experiment_name}_{str(N)}_params.txt",
        float,
    )
    n_params = len(params)

    #############################
    # Loop through all params
    #############################
    # Arrays to save results
    sims_mean = np.zeros(n_params)
    sims_sd = np.zeros(n_params)
    sims_range = np.zeros(n_params)
    if store_sims:
        sims_out = np.zeros((n_params, len(tas_in)))

    if eval == "SMAP":
        sims_rmse = np.zeros(n_params)
        sims_ubrmse = np.zeros(n_params)
    elif eval == "NLDAS":
        sims_rmse = [np.zeros(n_params) for _ in range(3)]
        sims_ubrmse = [np.zeros(n_params) for _ in range(3)]

    for iparam in range(n_params):
        # Read in correct order!
        Ts = params[iparam][0]
        Tm = params[iparam][1]
        wiltingp = params[iparam][2]
        awCap = params[iparam][3]
        alpha = params[iparam][4]
        betaHBV = params[iparam][5]
        if Kpet_name == "crop":
            GS_start = int(params[iparam][6])
            GS_length = int(params[iparam][7])
            L_ini = params[iparam][8]
            L_dev = params[iparam][9]
            L_mid = params[iparam][10]
            L_late = 1.0 - (L_ini + L_dev + L_mid)
            Kc_ini = params[iparam][11]
            Kc_mid = params[iparam][12]
            Kc_end = params[iparam][13]
            Kmin = params[iparam][14]
            Kmax = params[iparam][15]
            c_lai = params[iparam][16]
            # Construct Kc timeseries
            Kpet_in = construct_Kpet_crop(
                GS_start,
                GS_length,
                L_ini,
                L_dev,
                L_mid,
                L_late,
                Kc_ini,
                Kc_mid,
                Kc_end,
                Kmin,
                Kmax,
                c_lai,
                lai_in,
            )
        elif Kpet_name == "gen":
            Kmin = params[iparam][6]
            Kmax = params[iparam][7]
            c_lai = params[iparam][8]
            # Construct Kc timeseries
            Kpet_in = construct_Kpet_gen(Kmin, Kmax, c_lai, lai_in)

        # Initial conditions
        Ws_init = wiltingp + (awCap / 2.0)  # Initial soil moisture
        Wi_init = 0.0  # Canopy water storage
        Sp_init = 0.0  # Snowpack

        # Assume 1m root depth
        rootDepth = 1.0

        # Run it
        sim = wbm_jax(
            tas=tas_in + tas_delta,
            prcp=prcp_in * prcp_factor,
            Kpet=Kpet_in,
            Ws_init=Ws_init,
            Wi_init=Wi_init,
            Sp_init=Sp_init,
            lai=lai_in,
            phi=phi,
            params=(Ts, Tm, awCap, wiltingp, rootDepth, alpha, betaHBV),
        )

        # Store metrics
        if store_sims:
            sims_out[iparam, :] = sim

        sims_mean[iparam] = jnp.mean(sim)
        sims_sd[iparam] = jnp.std(sim)
        sims_range[iparam] = jnp.max(sim) - jnp.min(sim)

        sim_centered = sim - jnp.mean(sim)
        if eval == "SMAP":
            sims_rmse[iparam] = jnp.sqrt(jnp.mean((sim - obs) ** 2))
            sims_ubrmse[iparam] = jnp.sqrt(
                jnp.mean((sim_centered - obs_centered) ** 2)
            )
        elif eval == "NLDAS":
            for io in range(len(obs_list)):
                sims_rmse[io][iparam] = jnp.sqrt(
                    jnp.mean((sim - obs[io]) ** 2)
                )
                sims_ubrmse[io][iparam] = jnp.sqrt(
                    jnp.mean((sim_centered - obs_centered[io]) ** 2)
                )

    # Store simulation + stats if desired
    out = {}
    if store_sims:
        out["sims"] = sims_out
        out["mean"] = sims_mean
        out["sd"] = sims_sd
        out["range"] = sims_range
        if eval == "SMAP":
            out["SMAP_rmse"] = sims_rmse
            out["SMAP_ubrmse"] = sims_ubrmse
        elif eval == "NLDAS":
            for io in range(len(obs_list)):
                out[f"{obs_list[io]}_rmse"] = sims_rmse[io]
                out[f"{obs_list[io]}_ubrmse"] = sims_ubrmse[io]

        np.savez(
            f"{project_data_path}/WBM/SA/{save_name}_sim_stats.npz", **out
        )

    ####################
    # Calculate indices
    ####################
    # Problem spec
    with open(f"{project_data_path}/WBM/SA/{experiment_name}.json") as f:
        params_dict = json.load(f)

    param_names = list(params_dict.keys())
    sp = ProblemSpec(
        {
            "num_vars": n_params,
            "names": param_names,
            "bounds": [params_dict[param] for param in param_names],
        }
    ).set_samples(params)

    # Analyze all
    df_out_total = []
    df_out_second = []

    # Set up list
    if eval == "SMAP":
        metrics = [sims_mean, sims_sd, sims_range, sims_rmse, sims_ubrmse]
        metric_names = ["mean", "sd", "range", "rmse_SMAP", "ubrmse_SMAP"]
    elif eval == "NLDAS":
        metrics = (
            [sims_mean, sims_sd, sims_range]
            + [rmse for rmse in sims_rmse]
            + [ubrmse for ubrmse in sims_ubrmse]
        )
        metric_names = (
            ["mean", "sd", "range"]
            + [f"rmse_{obs}" for obs in obs_list]
            + [f"ubrmse_{obs}" for obs in obs_list]
        )
    else:
        metrics = [sims_mean, sims_sd, sims_range]
        metric_names = ["mean", "sd", "range"]

    # Loop and calculate
    for metric, metric_name in zip(metrics, metric_names):
        # Calculate
        sp.set_results(metric)
        sp.analyze_sobol()
        # Store
        total, first, second = sp.to_df()

        df_tmp = pd.merge(total, first, left_index=True, right_index=True)
        df_tmp["metric"] = metric_name
        df_out_total.append(
            df_tmp.reset_index().rename(columns={"index": "param"})
        )

        second = second.reset_index().rename(columns={"index": "param"})
        second["metric"] = metric_name
        second["param1"] = second["param"].apply(lambda x: x[0])
        second["param2"] = second["param"].apply(lambda x: x[1])
        df_out_second.append(second.drop(columns="param"))

    # Save
    df_out_total = pd.concat(df_out_total)
    df_out_second = pd.concat(df_out_second)

    df_out_total.to_csv(
        f"{project_data_path}/WBM/SA/{save_name}_res_total.csv", index=False
    )
    df_out_second.to_csv(
        f"{project_data_path}/WBM/SA/{save_name}_res_2order.csv", index=False
    )
