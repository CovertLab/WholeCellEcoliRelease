#! /usr/bin/env python

"""
NCA methods to use to solve E = AP given E and specified network connections in A.
"""

from typing import Callable, cast, Dict, List, Set, Tuple

import numpy as np
import scipy.linalg
import scipy.optimize
import swiglpk as glp

from wholecell.utils import parallelization


# Function names of the different NCA methods that have been implemented below
METHODS = ['robust_nca', 'constrained_nca', 'fast_nca', 'random_nca']


class NotValidMatrixError(Exception):
    pass


def nonnegative_least_squares(A: np.ndarray, B: np.ndarray) -> np.ndarray:
    """
    Solve nonnegative least squares with two matrices.
    min ||AX - B||
     st X >= 0
    """

    X = np.zeros((A.shape[1], B.shape[1]))
    for i, b in enumerate(B.T):
        X[:, i] = scipy.optimize.nnls(A, b)[0]

    return X

def nca_criteria_check(A: np.ndarray, tfs: np.ndarray, verbose: bool = True) -> Tuple[np.ndarray, np.ndarray]:
    """
    Check criteria for the A matrix of NCA (E = AP).
    - full column rank in A
    - removing column and associated rows leads to rank of m-1 for m columns in A

    Args:
        A: initial A matrix for NCA algorithm (n genes, m TFs)
        tfs: IDs of transcription factors for each column in A
        verbose: If set, prints updates about removed columns

    Returns:
        A: updated A matrix with columns removed to satisfy NCA criteria
        tfs: updated TF IDs corresponding to the new A
    """

    def nan_rank(M: np.ndarray) -> int:
        """
        Get the matrix rank for a matrix that can have NaN values by converting
        NaN to 1. NaN can be in the A matrix to identify uncertainty in the
        regulation direction so it should take on a real value when calculating
        rank.
        """

        mask = ~np.isfinite(M)
        if np.any(mask):
            N = M.copy()
            N[mask] = 1
        else:
            N = M

        try:
            rank = np.linalg.matrix_rank(N)
        except ValueError as e:
            raise NotValidMatrixError(f'Could not get the rank of the matrix: {e!r}')

        return rank

    if verbose:
        print(f'A shape: {A.shape}')

    # Filter out TFs that did not match any known genes for regulation (empty columns)
    col_sum = np.sum(A != 0, axis=0)
    if np.any(col_sum == 0):
        tf_idx_with_regulation = np.where(col_sum != 0)[0]
        if verbose:
            print('Removing empty columns for TFs:')
            for tf in tfs[col_sum == 0]:
                print(f'\t{tf}')
        A, tfs = nca_criteria_check(A[:, tf_idx_with_regulation], tfs[tf_idx_with_regulation], verbose=verbose)

    n_cols = A.shape[1]

    # Check if not full rank A matrix
    rank = nan_rank(A)
    if rank < n_cols:
        if verbose:
            print('A matrix is not full rank because of TFs:')
        min_entries = np.inf
        removed_tf = None
        mask = np.ones(n_cols, bool)
        for i in range(n_cols):
            col_mask = np.ones(n_cols, bool)
            col_mask[i] = False
            sub_rank = nan_rank(A[:, col_mask])
            if sub_rank == rank:
                n_entries = np.sum(A[:, i] != 0)
                if n_entries < min_entries:
                    mask = col_mask
                    min_entries = n_entries
                    removed_tf = tfs[i]
                if verbose:
                    print(f'\t{tfs[i]}')
        if verbose:
            print(f'\tRemoved {removed_tf}')

        A, tfs = nca_criteria_check(A[:, mask], tfs[mask], verbose=verbose)
        rank = nan_rank(A)
        n_cols = A.shape[1]

    # Check that reduced matrices have expected rank
    for i in range(n_cols):
        row_mask = A[:, i] == 0
        col_mask = np.ones(n_cols, bool)
        col_mask[i] = False

        reduced_A = A[row_mask, :][:, col_mask]
        sub_rank = nan_rank(reduced_A)
        if sub_rank != rank - 1:
            # Filter out empty columns with rows removed
            col_sum = np.sum(A[row_mask, :] != 0, axis=0)
            tf_idx_with_regulation = col_sum != 0
            tf_idx_with_regulation[i] = True
            if not np.all(tf_idx_with_regulation):
                # TODO: turn into function with same code above
                if verbose:
                    print(f'Removing empty columns from reduced matrix of {tfs[i]} for TFs:')
                    for tf in tfs[~tf_idx_with_regulation]:
                        print(f'\t{tf}')
                A, tfs = nca_criteria_check(A[:, tf_idx_with_regulation], tfs[tf_idx_with_regulation], verbose=verbose)
            else:
                # TODO; turn into function with same code above
                if verbose:
                    print('Reduced matrix is not full rank because of TFs:')
                rows_removed = A[row_mask, :]
                min_entries = np.inf
                removed_tf = None
                mask = np.ones(n_cols, bool)
                for j in range(n_cols):
                    col_mask = np.ones(n_cols, bool)
                    col_mask[j] = False
                    new_rank = nan_rank(rows_removed[:, col_mask])
                    if new_rank == sub_rank:
                        n_entries = np.sum(A[:, j] != 0)
                        if n_entries < min_entries:
                            mask = col_mask
                            min_entries = n_entries
                            removed_tf = tfs[j]
                        if verbose:
                            print(f'\t{tfs[j]}')
                if verbose:
                    print(f'\tRemoved {removed_tf}')

                A, tfs = nca_criteria_check(A[:, mask], tfs[mask], verbose=verbose)
            break

    return A, tfs

def random_nca(E: np.ndarray, A: np.ndarray, verbose: bool = True, **options) -> Tuple[np.ndarray, np.ndarray]:
    """
    Randomly assign values to A and solve for P with least squares to
    determine performance of randomly assigned A.

    Args:
        E: data to solve NCA for (n genes, m conditions)
        A: network connectivity (n genes, o TFs)
        verbose: if True, print progress updates
        options: solver options that are implemented for other methods

    Returns:
        A_est: estimated A based fit to data (n genes, o TFs)
        P_est: estimated P based fit to data (o TFs, m conditions)
    """

    if verbose:
        print('Solving with random A values...')

    # Set entries between -1.1 and -0.9 or 0.9 and 1.1 depending on sign of A
    A_est = cast(np.ndarray, np.sign(A))
    nonzero_mask = A_est != 0
    n_nonzero = np.sum(nonzero_mask)
    A_est[nonzero_mask] += 0.1 * (np.random.rand(n_nonzero) - 0.5)

    # Set ambiguous regulation between -1 and 1
    not_finite_mask = ~np.isfinite(A_est)
    n_not_finite = np.sum(not_finite_mask)
    A_est[not_finite_mask] = np.random.rand(n_not_finite) * 2 - 1

    # Solve for P
    P_est = np.linalg.lstsq(A_est, E, rcond=None)[0]

    return A_est, P_est


def fast_nca(E: np.ndarray,
        A: np.ndarray,
        verbose: bool = True,
        status_step: float = 0.1,
        **options) -> Tuple[np.ndarray, np.ndarray]:
    """
    Perform FastNCA on dataset E with network connectivity specified by A for the
    problem: E = AP. Based on matlab implementation from Chang et al. 2008.
    https://www.eee.hku.hk/~cqchang/FastNCA.m

    Args:
        E: data to solve NCA for (n genes, m conditions)
        A: network connectivity (n genes, o TFs)
        verbose: if True, print progress updates
        status_step: print status update every time this fraction of steps has
            been completed
        options: solver options that are implemented for other methods

    Returns:
        A_est: estimated A based fit to data (n genes, o TFs)
        P_est: estimated P based fit to data (o TFs, m conditions)
    """

    if verbose:
        print('Solving with FastNCA...')
    n_cols = A.shape[1]
    U, S, V = scipy.linalg.svd(E)
    U = U[:, :n_cols]
    S = np.diag(S[:n_cols])
    V = V[:n_cols, :]
    E_approx = U.dot(S).dot(V)

    A_est = np.zeros_like(A)
    next_update = status_step
    for i in range(n_cols):
        if verbose and (i + 1) / n_cols >= next_update:
            complete = np.floor((i + 1) / n_cols * 100)
            print(f'{complete:.0f}% complete...')
            next_update += status_step

        U0 = U[A[:, i] != 0, :]
        M, _, _ = scipy.linalg.svd(U0)
        v = M[:, 0]

        A_est[A[:, i] != 0, i] = v / np.sum(np.abs(v)) * np.sum(A[:, i] != 0)

    P_est = np.linalg.lstsq(A_est, E_approx, rcond=None)[0]

    return A_est, P_est

def robust_nca(
        E: np.ndarray,
        A: np.ndarray,
        verbose: bool = True,
        status_step: float = 0.1,
        n_iters: int = 5,
        error_tolerance = 1e-3,
        **options) -> Tuple[np.ndarray, np.ndarray]:
    """
    Perform ROBNCA on dataset E with network connectivity specified by A for the
    problem: E = AP. Based on method in Noor et al. Bioinformatics. 2013.

    Args:
        E: data to solve NCA for (n genes, m conditions)
        A: network connectivity (n genes, o TFs)
        verbose: if True, print progress updates
        status_step: print status update every time this fraction of steps has
            been completed
        n_iters: maximum number of iterations before exiting
        error_tolerance: difference in error for early stopping
        options: solver options that are implemented for other methods

    Returns:
        A_est: estimated A based fit to data (n genes, o TFs)
        P_est: estimated P based fit to data (o TFs, m conditions)
    """

    if verbose:
        print('Solving with ROBNCA...')
    n_genes = E.shape[0]
    n_tfs = A.shape[1]
    A_est = A.astype(float)
    outliers = np.zeros_like(E)
    zero_mask = A == 0
    A_est[~np.isfinite(A_est)] = 1

    lambda_ = 2  # parameter that can be tuned for sparsity to determine outliers
    next_update = status_step
    old_error = np.inf
    for it in range(n_iters):
        # Update S
        X = E - outliers
        S = np.linalg.inv(A_est.T.dot(A_est)).dot(A_est.T).dot(X)

        # Update A
        Qinv = np.linalg.inv(S.dot(S.T))
        for col, a in enumerate(A_est):
            zero_idx = np.where(a == 0)[0]
            Cn = np.zeros((len(zero_idx), n_tfs))
            Cn[range(len(zero_idx)), zero_idx] = 1
            psi = Cn.dot(Qinv).dot(Cn.T)
            w = S.dot(X.T[:, col])
            A_est[col, :] = Qinv.dot(w - Cn.T.dot(np.linalg.inv(psi)).dot(Cn).dot(Qinv).dot(w))
        A_est[zero_mask] = 0

        # Update outliers
        E_est = A_est.dot(S)
        B = E - E_est
        norm = np.linalg.norm(B, axis=0) / np.sqrt(n_genes)
        outliers = B * np.fmax(0, (norm - lambda_ / 2)) / norm

        # Calculate error for early break
        mask = E != 0  # filter to remove inf and nan when expression has been adjusted to 0 in ISNCA
        error = np.sqrt(np.mean((B[mask] / E[mask])**2))
        if (old_error - error) / error < error_tolerance:
            print(f'Completed after {it+1} iterations')
            break

        # Print progress update for certain iterations
        if verbose and (it + 1) / n_iters >= next_update:
            complete = np.floor((it + 1) / n_iters * 100)
            print(f'{complete:.0f}% complete (error = {error:.4f})...')
            next_update += status_step

    P_est = S
    return A_est, P_est

def constrained_nca(E: np.ndarray, A: np.ndarray, verbose: bool = True, **options) -> Tuple[np.ndarray, np.ndarray]:
    """
    Perform constrained NCA on dataset E with network connectivity specified by
    A for the problem: E = AP. Based on method in Chang et al. ICA. 2009.
    Original method was only for non-negative constraint but this has been
    modified to accept non-negative and non-positive constraints.

    Args:
        E: data to solve NCA for (n genes, m conditions)
        A: network connectivity (n genes, o TFs), sign should indicate
            sign constraint for problem solution (eg. positive entries will
            stay positive, negative entries will stay negative)
        verbose: if True, print progress updates
        options: solver options that are implemented for other methods

    Returns:
        A_est: estimated A based fit to data (n genes, o TFs)
        P_est: estimated P based fit to data (o TFs, m conditions)
    """

    def int_array(array):
        ia = glp.intArray(len(array) + 1)
        ia[0] = -1
        for (i, value) in enumerate(array):
            ia[i + 1] = int(value)
        return ia

    def double_array(array):
        da = glp.doubleArray(len(array) + 1)
        da[0] = np.nan
        for (i, value) in enumerate(array):
            da[i + 1] = float(value)
        return da

    if verbose:
        print('Solving constrained NCA...')
    A_cols = A.shape[1]
    U, _, _ = scipy.linalg.svd(E)
    C = U[:, A_cols:]
    C_cols = C.shape[1]

    # Create LP objects
    lp = glp.glp_create_prob()
    smcp = glp.glp_smcp()
    glp.glp_init_smcp(smcp)
    if not verbose:
        smcp.msg_lev = glp.GLP_MSG_ERR

    # Set up variables
    t_idx = 1  # +1 offset for GLPK
    n_entries = int(np.sum(A != 0))
    glp_rows = A_cols * (2*C_cols + 1)
    glp.glp_add_rows(lp, glp_rows)
    glp_cols = n_entries + 1
    glp.glp_add_cols(lp, glp_cols)

    # Add constraints
    glp.glp_set_col_bnds(lp, t_idx, glp.GLP_FR, 0, 0)
    bound = 0.
    row_idx = 1  # GLPK +1 offset
    col_idx = 2  # Start at 2 for GLPK +1 offset and t_idx
    result_mapping = {}
    for col, data in enumerate(A.T):
        nonzero_idx = np.where(data)[0]
        n_nonzero = len(nonzero_idx)
        data_sign = np.sign(data[nonzero_idx])

        ## -t < ca < t
        for i in range(C_cols):
            values = np.hstack((1, C[nonzero_idx, i]))
            length = len(values)
            col_idxs = int_array(np.hstack((t_idx, range(col_idx, col_idx + length))))

            # ca + t > 0 (-t < ca)
            glp.glp_set_mat_row(lp, row_idx, length, col_idxs, double_array(values))
            glp.glp_set_row_bnds(lp, row_idx, glp.GLP_LO, bound, bound)
            row_idx += 1

            # ca - t < 0 (ca < t)
            values[0] = -1  # update t sign
            glp.glp_set_mat_row(lp, row_idx, length, col_idxs, double_array(values))
            glp.glp_set_row_bnds(lp, row_idx, glp.GLP_UP, bound, bound)
            row_idx += 1

        ## sum(a) = Lj
        col_idxs = int_array(range(col_idx, col_idx + n_nonzero))
        values = double_array(data_sign)
        glp.glp_set_mat_row(lp, row_idx, n_nonzero, col_idxs, values)
        glp.glp_set_row_bnds(lp, row_idx, glp.GLP_FX, float(n_nonzero), float(n_nonzero))
        row_idx += 1

        ## A sign constraint
        for i, sign in enumerate(data_sign):
            glp_col = col_idx + i
            result_mapping[glp_col - 1] = (nonzero_idx[i], col)
            if not np.isfinite(sign):
                bound_type = glp.GLP_FR
            elif sign > 0:
                bound_type = glp.GLP_LO
            else:
                bound_type = glp.GLP_UP
            glp.glp_set_col_bnds(lp, glp_col, bound_type, bound, bound)
        col_idx += n_nonzero

    # Set objective
    glp.glp_set_obj_coef(lp, t_idx, 1)

    # Solve LP
    glp.glp_set_obj_dir(lp, glp.GLP_MIN)
    result = glp.glp_simplex(lp, smcp)
    if result != 0:
        raise RuntimeError('Could not solve LP problem')
    solution = glp.get_col_primals(lp)

    # Recreate A matrix from LP solution
    A_est = np.zeros_like(A)
    for idx, (i, j) in result_mapping.items():
        A_est[i, j] = solution[idx]

    # Solve for P
    P_est = np.linalg.lstsq(A_est, E, rcond=None)[0]

    return A_est, P_est

def iterative_sub_nca(
        method: Callable[[np.ndarray, np.ndarray], Tuple[np.ndarray, np.ndarray]],
        E: np.ndarray,
        A: np.ndarray,
        tfs: np.ndarray,
        statistics: Callable[[np.ndarray, np.ndarray, np.ndarray, np.ndarray, Dict[str, Dict[str, int]], np.ndarray], None] = None,
        statistics_args: Tuple = (),
        n_iters: int = 100,
        splits: int = 2,
        parallel: bool = True,
        robust_iters: int = 1,
        status_step: float = 0.1,
        verbose: bool = False,
        ) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Iterative sub-network component analysis method applied to any NCA method
    above. Based on method in Jayavelu et al. BMC Bioinformatics. 2015.

    Args:
        method: NCA method implemented in this file that takes E and A as args
        E: data to solve NCA for (n genes, m conditions)
        A: network connectivity (n genes, o TFs)
        tfs: IDs of transcription factors for each column in A
        statistics: function to provide match statistics at each iteration,
            needs to take E, A, P, tfs as the first args
        statistics_args: args to the statistics function that are not determined
            in this method
        n_iters: maximum number of iterations to perform
        splits: maximum number of sub-networks to split into
        parallel: if True, runs each sub-network NCA in parallel
        robust_iters: maximum number of iterations for robust_nca
        status_step: print status update every time this fraction of steps has
            been completed
        verbose: if True, displays additional information

    Returns:
        A_est: estimated A based fit to data (n genes, o TFs)
        P_est: estimated P based fit to data (o TFs, m conditions)
    """

    def divide_network(
            E: np.ndarray,
            A: np.ndarray,
            tfs: np.ndarray,
            verbose: bool,
            max_divisions: int,
            ) -> Tuple[List[np.ndarray], List[np.ndarray], List[np.ndarray], List[Set[int]], List[Set[int]]]:
        """
        Divide the network into subnetworks with unique and common TFs (eq. 3 and 4).

        Args:
            E: original expression matrix (n genes, o TFs)
            A: original network topology mapping (o TFs, m conditions)
            tfs: original TF IDs corresponding to columns of A
            verbose: If set, prints updates about removed columns
            max_divisions: the maximum number of times to subdivide the A matrix

        Returns:
            E_divided: reduced E matrices associated with each subnetwork,
                only include genes produced by corresponding A mapping
            A_divided: reduced A matrices assocaited with each subnetwork,
                only include genes and TFs with relationships that satisfy
                the NCA criteria
            tfs_divided: TF IDs associated withe the columns of each reduced
                A matrix for each subnetwork
            common_genes: gene indices (rows in original E and A matrices)
                shared with this subnetwork and any others
            unique_genes: gene indices (rows in original E and A matrices)
                unique to each subnetwork and not shared with any others
        """

        E_divided = []  # type: List[np.ndarray]
        A_divided = []
        tfs_divided = []
        divided_genes = []
        common_genes = []
        unique_genes = []

        n_divisions = 0
        removed_tfs = set(tfs)
        while removed_tfs and n_divisions < max_divisions:
            n_divisions += 1
            reduced_A = np.array([tf in removed_tfs for tf in tfs])
            try:
                Ai, tfsi = nca_criteria_check(A[:, reduced_A], tfs[reduced_A], verbose=verbose)
            except NotValidMatrixError as e:
                print(f'Warning: could only make {len(E_divided)} divisions: {e!r}')
                break
            removed_tfs = removed_tfs - set(tfsi)

            regulated_genes = np.unique(np.where(Ai)[0])
            E_divided.append(E[regulated_genes, :])
            A_divided.append(Ai[regulated_genes, :])
            tfs_divided.append(tfsi)
            divided_genes.append(set(regulated_genes))

        for i, genes in enumerate(divided_genes):
            others = {g for j, genes in enumerate(divided_genes) if i != j for g in genes}
            common_genes.append(genes.intersection(others))
            unique_genes.append(genes.difference(others))

        return E_divided, A_divided, tfs_divided, common_genes, unique_genes

    def solve_networks(
            method: Callable,
            E_divided: List[np.ndarray],
            A_divided: List[np.ndarray],
            cpus: int,
            **kwargs) -> Tuple[List[np.ndarray], List[np.ndarray]]:
        """Solve the divided networks with the given method."""

        A_hat = []
        P_hat = []

        cpus = min(cpus, len(E_divided))
        if cpus > 1:
            pool = parallelization.pool(cpus)
            results = [
                pool.apply_async(method, (E, A), kwds=kwargs)
                for E, A in zip(E_divided, A_divided)
                ]
            pool.close()
            pool.join()

            for result in results:
                A_est, P_est = result.get()
                A_hat.append(A_est)
                P_hat.append(P_est)
        else:
            for E, A in zip(E_divided, A_divided):
                A_est, P_est = method(E, A, **kwargs)
                A_hat.append(A_est)
                P_hat.append(P_est)

        return A_hat, P_hat

    def assemble_network(
            n_genes: int,
            n_tfs: int,
            A_hat: List[np.ndarray],
            P_hat: List[np.ndarray],
            common_genes: List[Set[int]],
            unique_genes: List[Set[int]],
            ) -> Tuple[np.ndarray, np.ndarray]:
        """
        Assemble subnetwork solutions into combined A and P solutions (eq. 8).

        Args:
            n_genes: number of genes in the reconstructed A matrix
            n_tfs: number of TFs in the reconstructed A matrix
            A_hat: estimated A for each subnetwork
            P_hat: estimated P for each subnetwork
            common_genes: indices of common genes shared with other subnetworks
                for each subnetwork
            unique_genes: indices of unique genes to each subnetwork

        Returns:
            A_est: assembled A matrix to estimate all subnetworks
            P_est: assembled P matrix to estimate all subnetworks
        """

        A_est = np.zeros((n_genes, n_tfs))
        last_col = 0
        for sub_A, common, unique in zip(A_hat, common_genes, unique_genes):
            gene_idx = np.sort(np.hstack((np.array(list(common)), np.array(list(unique))))).astype(int)
            start_col = last_col
            last_col = start_col + sub_A.shape[1]
            A_est[gene_idx, start_col:last_col] = sub_A

        P_est = np.vstack(P_hat)

        return A_est, P_est

    def update_E(
            E: np.ndarray,
            A_hat: List[np.ndarray],
            P_hat: List[np.ndarray],
            common_genes: List[Set[int]],
            unique_genes: List[Set[int]],
            attenuation: float,
            ) -> List[np.ndarray]:
        """
        Update E submatrices for the next iteration (eq. 11).

        Args:
            E: original expression matrix
            A_hat: estimated A for each subnetwork
            P_hat: estimated P for each subnetwork
            common_genes: indices of common genes shared with other subnetworks
                for each subnetwork
            unique_genes: indices of unique genes to each subnetwork
            attenuation: fraction to multiply shared predictions by to reduce
                E for each subnetwork

        Returns:
            E_divided: updated E matrices for each subnetwork by subtracting out
                contributions from other subnetworks to shared components
        """

        # Contributions from each subnetwork
        Ts = []
        for A, P in zip(A_hat, P_hat):
            Ts.append(A.dot(P))

        # Update each subnetwork's E matrix with contributions from other subnetworks
        E_divided = []
        for i, (common, unique) in enumerate(zip(common_genes, unique_genes)):
            gene_idx = np.sort(np.hstack((np.array(list(common)), np.array(list(unique))))).astype(int)
            new_E = E.copy()[gene_idx, :]

            # Map original index (global among all subnetworks) to reduced index (row in E)
            update_idx = {original_idx: reduced_idx
                for reduced_idx, original_idx in enumerate(gene_idx)
                if original_idx in common}

            # Iterate through every other subnetwork to update common contributions
            for j, (T, T_common, T_unique) in enumerate(zip(Ts, common_genes, unique_genes)):
                if i == j:
                    continue

                T_idx = np.sort(np.hstack((np.array(list(T_common)), np.array(list(T_unique))))).astype(int)
                for row, idx in zip(T, T_idx):
                    if idx in common:
                        # TODO: handle values <0?
                        new_E[update_idx[idx], :] -= attenuation * row

            # Prevent negative values for better stability since this is not biologically possible
            if np.any(new_E < 0):
                print('Values in new_E adjusted to 0 from negative values.')
                new_E[new_E < 0] = 0

            E_divided.append(new_E)

        return E_divided

    # Parameters for approach
    error_threshold = 1e-5
    old_error = np.inf
    best_error = np.inf
    best_A = np.array([])
    best_P = np.array([])

    E_divided, A_divided, tfs_divided, common_genes, unique_genes = divide_network(E, A, tfs, verbose, splits)
    new_tfs = np.array([tf for tfs in tfs_divided for tf in tfs])
    n_genes = A.shape[0]
    n_tfs = len(new_tfs)

    for it in range(n_iters):
        # Solve for A and P in subnetworks and overall problem
        try:
            A_hat, P_hat = solve_networks(method, E_divided, A_divided, splits if parallel else 1,
                n_iters=robust_iters, status_step=status_step, verbose=verbose)
            A_est, P_est = assemble_network(n_genes, n_tfs, A_hat, P_hat, common_genes, unique_genes)
        except Exception as e:
            print(f'{type(e).__name__} while running ISNCA: {e.args[0]}')
            break
        except KeyboardInterrupt:
            break

        # Check if solution has converged
        E_est = A_est.dot(P_est)
        predicted_mask = E_est.sum(1) != 0
        error = np.linalg.norm(E[predicted_mask, :] - E_est[predicted_mask, :])
        # TODO: better criteria for exiting
        if np.abs(error - old_error) / error < error_threshold:
            break
        old_error = error
        print(f'Iteration {it} error: {error:.3f}')

        # Print statistics for this iteration
        if verbose and statistics:
            statistics(E, A_est, P_est, new_tfs, *statistics_args)

        # Check if this is the best solution so far and store results
        if error < best_error:
            best_error = error
            best_A = A_est.copy()
            best_P = P_est.copy()

        # Update E_divided matrices
        attenuation = min((it + 1) / (10*len(E_divided)), 1)  # between 0 and 1
        E_divided = update_E(E, A_hat, P_hat, common_genes, unique_genes, attenuation)

    print(f'Completed ISNCA with low error of {best_error:.3f}')

    return best_A, best_P, new_tfs
