from datetime import datetime

import jax
import jax.numpy as jnp
import numpy as np
import optax

from prediction import make_prediction
from read_inputs import read_hindcast_inputs
from utils.global_paths import project_data_path
from utils.initial_params import (
    constants,
    initial_params,
)


def test():
    print(project_data_path)


def train_and_store(
    subset_name,
    obs_name,
    _error_fn,
    error_fn_name,
    initial_theta,
    params_lower,
    params_upper,
    param_names,
    n_epochs_min=10,
    n_epochs_max=30,
    patience=5,
    batch_size=2**5,
    opt="adam",
    learning_rate=1e-3,
    val_inds=[],
    reg_const=0.001,
):
    #############################################
    # Loss functions
    ############################################
    # Prediction loss
    def prediction_loss(
        theta, constants, x_forcing_nt, x_forcing_nyrs, x_maps, ys
    ):
        prediction = make_prediction(
            theta, constants, x_forcing_nt, x_forcing_nyrs, x_maps
        )
        return _error_fn(prediction, ys)

    # Regularization loss
    def reg_loss(theta, initial_params, params_lower, params_upper):
        return jnp.nansum(
            (theta - initial_params) ** 2
            / ((theta - params_lower) * (params_upper - theta))
        )

    # Total loss
    def loss_fn(
        theta,
        reg_const,
        initial_params,
        params_lower,
        params_upper,
        constants,
        x_forcing_nt,
        x_forcing_nyrs,
        x_maps,
        ys,
    ):
        return prediction_loss(
            theta, constants, x_forcing_nt, x_forcing_nyrs, x_maps, ys
        ) + reg_const * reg_loss(
            theta, initial_params, params_lower, params_upper
        )

    # jit and vmap it
    pred_loss_value = jax.jit(
        jax.vmap(prediction_loss, in_axes=(None, None, 0, 0, 0, 0), out_axes=0)
    )
    loss_value_and_grad = jax.jit(
        jax.vmap(
            jax.value_and_grad(loss_fn),
            in_axes=(None, None, None, None, None, None, 0, 0, 0, 0),
            out_axes=0,
        )
    )

    ###########################
    # Setup
    ###########################
    # Read data
    ys, x_forcing_nt, x_forcing_nyrs, x_maps = read_hindcast_inputs(
        subset_name, obs_name, True
    )
    N = ys.shape[0]

    # Get train/val split over space
    if len(val_inds) > 0:
        ys_val, x_forcing_nt_val, x_forcing_nyrs_val, x_maps_val = (
            ys[val_inds],
            x_forcing_nt[val_inds],
            x_forcing_nyrs[val_inds],
            x_maps[val_inds],
        )
        train_inds = np.array([n for n in np.arange(N) if n not in val_inds])
        ys_train, x_forcing_nt_train, x_forcing_nyrs_train, x_maps_train = (
            ys[train_inds],
            x_forcing_nt[train_inds],
            x_forcing_nyrs[train_inds],
            x_maps[train_inds],
        )
    else:
        ys_train, x_forcing_nt_train, x_forcing_nyrs_train, x_maps_train = (
            ys,
            x_forcing_nt,
            x_forcing_nyrs,
            x_maps,
        )

    # Memory management
    del ys, x_forcing_nt, x_forcing_nyrs, x_maps

    # Define mini-batch hyper-parameters
    N_train = ys_train.shape[0]
    n_minibatches = int(np.ceil(N_train / batch_size))

    # Initial parameters
    theta = initial_theta

    # Optimizer
    if opt == "adam":
        adam = optax.adam(learning_rate=learning_rate)
        opt_fn = adam.update
        opt_state = adam.init(theta)
    elif opt == "sgd":
        learning_rate = 1e-5
        opt_state = None

        def sgd(gradients, state):
            return -learning_rate * gradients, state

        opt_fn = sgd

    # Loss
    train_loss_out = np.empty(n_epochs_max + 1)
    pred_loss_out = np.empty(n_epochs_max + 1)
    reg_loss_out = np.empty(n_epochs_max + 1)
    val_loss_out = np.empty(n_epochs_max + 1)

    # Early stopping
    best_loss = jnp.inf

    #####################
    # Initial results
    #####################
    datetime_str = datetime.now().strftime("%Y%m%d-%H%M")

    # Used to discern different starting values
    random_str = str(abs(theta[0])).replace(".", "")[:5]
    training_name = f"{error_fn_name}_{random_str}r"
    if len(val_inds) > 0:
        training_name += f"_{str(val_inds[0])}val"

    # Write initial results
    out_file_path = f"{project_data_path}/WBM/calibration/{subset_name}/{obs_name}/training_res/{training_name}_{datetime_str}.txt"
    f = open(out_file_path, "w")
    f.write(
        f"epoch metric train_loss pred_loss reg_loss val_loss {' '.join(param_names)}\n"
    )

    pred_loss_init = jnp.mean(
        pred_loss_value(
            theta,
            constants,
            x_forcing_nt_train,
            x_forcing_nyrs_train,
            x_maps_train,
            ys_train,
        )
    )
    if len(val_inds) > 0:
        val_loss_init = jnp.mean(
            pred_loss_value(
                theta,
                constants,
                x_forcing_nt_val,
                x_forcing_nyrs_val,
                x_maps_val,
                ys_val,
            )
        )
    else:
        val_loss_init = -1.0

    reg_loss_init = reg_loss(theta, initial_params, params_lower, params_upper)
    train_loss_init = pred_loss_init + (reg_const * reg_loss_init)

    print(
        f"Epoch 0 train loss: {train_loss_init:.4f} pred loss: {pred_loss_init:.4f}, reg_loss: {reg_loss_init:.4f}, val loss: {val_loss_init:.4f}"
    )
    theta_str = [str(param) for param in theta]
    f.write(
        f"0 {error_fn_name} {train_loss_init:.4f} {pred_loss_init:.4f} {reg_loss_init:.4f} {val_loss_init:.4f} {' '.join(theta_str)}\n"
    )

    ###########################
    # Training loop
    ###########################
    invalid_theta_count = 0
    stagnation_counter = 0

    for epoch in range(n_epochs_max):
        # Shuffle indices
        shuffled_inds = np.random.permutation(N_train)

        # Generate a mini-batch
        minibatch_inds = [
            shuffled_inds[(i * batch_size) : ((i + 1) * batch_size)]
            for i in range(n_minibatches)
        ]

        # For batch loss
        batch_loss = [None] * n_minibatches

        for idx, inds in enumerate(minibatch_inds):
            # Calculate gradient of loss function, update parameters
            loss, grads = loss_value_and_grad(
                theta,
                reg_const,
                initial_params,
                params_lower,
                params_upper,
                constants,
                x_forcing_nt_train[inds],
                x_forcing_nyrs_train[inds],
                x_maps_train[inds],
                ys_train[inds],
            )
            updates, opt_state = opt_fn(jnp.nanmean(grads, axis=0), opt_state)
            theta = optax.apply_updates(theta, updates)
            batch_loss[idx] = loss
            # Break if theta steps outside bounds
            if (theta < params_lower).any() or (theta > params_upper).any():
                # print("Found invalid parameter... re-initializaing")
                f.write("Invalid-parameter\n")
                theta = initial_theta
                invalid_theta_count += 1
                if invalid_theta_count > 2:
                    f.close()
                    return (
                        train_loss_out,
                        pred_loss_out,
                        reg_loss_out,
                        val_loss_out,
                        theta,
                    )
                else:
                    theta = np.random.uniform(
                        low=params_lower, high=params_upper
                    )
            # Break if theta goes to NaN
            if jnp.isnan(theta).any():
                # print("Found NaN parameter")
                f.write("NaN-parameter")
                f.close()
                return (
                    train_loss_out,
                    pred_loss_out,
                    reg_loss_out,
                    val_loss_out,
                    theta,
                )

        # Save all losses
        train_loss = jnp.nanmean(
            jnp.array([item for row in batch_loss for item in row])
        )
        train_loss_out[epoch] = train_loss
        reg_loss_out[epoch] = reg_loss(
            theta, initial_params, params_lower, params_upper
        )
        pred_loss = train_loss_out[epoch] - (reg_const * reg_loss_out[epoch])
        pred_loss_out[epoch] = pred_loss

        ###### Early stopping criteria
        if len(val_inds) > 0:
            val_loss = jnp.mean(
                pred_loss_value(
                    theta,
                    constants,
                    x_forcing_nt_val,
                    x_forcing_nyrs_val,
                    x_maps_val,
                    ys_val,
                )
            )
            val_loss_out[epoch] = val_loss
            # Early stopping for val loss based on binary
            if val_loss < best_loss:
                best_loss = val_loss
                stagnation_counter = 0
            else:
                stagnation_counter += 1
        else:
            val_loss_out[epoch] = -1.0
            # Early stopping for training loss based on abs pct diff
            if epoch > 0:
                pred_loss_prev = pred_loss_out[epoch - 1]
                pct_diff = jnp.abs(pred_loss - pred_loss_prev) / pred_loss_prev
                if (pred_loss < pred_loss_prev) and (pct_diff < 0.01):
                    stagnation_counter += 1
                else:
                    stagnation_counter = 0

        # Write every epoch
        theta_str = [str(param) for param in theta]
        f.write(
            f"{str(epoch + 1)} {error_fn_name} {train_loss_out[epoch]:.4f} {pred_loss_out[epoch]:.4f} {reg_loss_out[epoch]:.4f} {val_loss_out[epoch]:.4f} {' '.join(theta_str)}\n"
        )

        # Print
        print(
            f"Epoch {str(epoch + 1)} total loss: {train_loss_out[epoch]:.4f}, pred loss: {pred_loss_out[epoch]:.4f}, reg_loss: {reg_loss_out[epoch]:.4f}, val loss: {val_loss_out[epoch]:.4f}"
        )

        # Early stopping
        if (
            (stagnation_counter > patience)
            & (epoch > n_epochs_min)
            & (epoch < n_epochs_max)
        ):
            # print("Early stopping!")
            f.close()
            return (
                train_loss_out,
                pred_loss_out,
                reg_loss_out,
                val_loss_out,
                theta,
            )

    f.close()
    return train_loss_out, pred_loss_out, reg_loss_out, val_loss_out, theta
