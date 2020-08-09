# -*- coding: utf-8 -*-
"""
    @author: Derek Macklin
    @organization: Covert Lab, Department of Bioengineering, Stanford University
    @date: Created 3/22/2013
    @LastEditors: Hwrn
    @LastEditTime: 2020-08-09 16:58:15
    @FilePath: /WholeCellEcoliRelease/wholecell/utils/linear_programming.py
    @Description:
        保证可替换性 | Wrapper over cvxopt to call glpk (can be extended in the future to use other solvers)
    @TODO:
"""

import cvxopt
import cvxopt.solvers
import numpy as np

def linearProgramming(
        maximizeFlag: str, f: np.ndarray,
        A: np.ndarray, b: np.ndarray,
        lb: np.ndarray =None, ub: np.ndarray =None,
        constraintTypes=["S"], variableTypes=["C"],
        options=None):
    """
        @description:
            解决线性规划问题:
                maximize f'(x)
                subject to
                    Ax == b | x 的等式约束
                    lb <= x <= ub | x 的不等式约束
        @example:
            (trivial solution)
            >>> linearProgramming(
            ...     "maximize", np.ones(1), np.ones((1,1)), np.ones(1), np.ones(1), np.ones(1), "S", "C", None)
            (array([[1.]]), {'status': 'optimal', 'x': <1x1 matrix, tc='d'>, 's': <2x1 matrix, tc='d'>, 'y': <1x1 matrix, tc='d'>, 'z': <2x1 matrix, tc='d'>, 'primal objective': -1.0, 'dual objective': -1.0, 'gap': 0.0, 'relative gap': 0.0, 'primal infeasibility': 0.0, 'dual infeasibility': 0.0, 'primal slack': 0.0, 'dual slack': -0.0, 'residual as primal infeasibility certificate': None, 'residual as dual infeasibility certificate': None})

        @TODO constrainTypes 和 variableTypes 是必填的, 调整到 lb 和 ub 前面?

        @param maximizeFlag 若需最大化, 参数设为 "maximize", 否则进行最小化
        @param f
        @param constraintTypes,
        @param variableTypes,
        @param options
        @return {type}
    """
    # Input validation
    assert type(A) == np.ndarray and len(A.shape) == 2, "Matrix A must be a non-empty numpy array."
    m, n = A.shape

    assert type(b) == np.ndarray and b.size == m, "Vector b must be a numpy array with same number of rows as A."
    if lb:
        assert type(lb) == np.ndarray and lb.size == n, "Vector lb must be a numpy array with dimension equal to the number of columns of A."
    else:
        lb = - np.ones(A.shape[1]) * np.Inf
    if ub:
        assert type(ub) == np.ndarray and ub.size == n, "Vector ub must be a numpy array with dimension equal to the number of columns of A."
    else:
        ub = np.ones(A.shape[1]) * np.Inf
    if not np.all(np.array(constraintTypes) == "S"):
        assert False, "At least one of your constraintTypes is not currently supported."
    if not np.all(np.array(variableTypes) == "C"):
        assert False, "At least one of your variableTypes is not currently supported."
    if options != None and options["glpk"] != None: # 实际上作者还没用到这个参数
        raise Exception("Passing GLPK options is not currently supported.")

    if maximizeFlag.lower() == "maximize":
        f = f * (-1.0)            # 否则变成指针 DO *****NOT***** DO "*=" AS THAT WILL PERFORM MODIFICATION IN PLACE

    G = np.concatenate([np.identity(n), -1 * np.identity(n)], axis = 0)
    h = np.concatenate([ub, -1.0 * lb], axis = 0)

    A_m = cvxopt.matrix(A)
    f_m = cvxopt.matrix(f)
    G_m = cvxopt.matrix(G)
    h_m = cvxopt.matrix(h)
    b_m = cvxopt.matrix(b)

    # Turn off screen output (http://abel.ee.ucla.edu/cvxopt/userguide/coneprog.html?highlight=solvers.lp#s-external)
    cvxopt.solvers.options["LPX_K_MSGLEV"] = 0

    sol = cvxopt.solvers.lp(f_m, G_m, h_m, A = A_m, b = b_m, solver = "glpk")

    return np.array(sol["x"]), sol


if __name__ == "__main__":
    import doctest
    doctest.testmod()

    print("if nothing happened, then OK")
