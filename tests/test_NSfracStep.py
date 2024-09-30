import math
import re
import subprocess

import pytest

number = "([0-9]+.[0-9]+e[+-][0-9]+)"


@pytest.mark.parametrize("solver", ["IPCS_ABCN", "BDFPC_Fast"])
@pytest.mark.parametrize("num_processors", [1])
def test_spatial_rate_of_convergence(num_processors, solver):
    # Simulation parameters
    p_err = []
    u0_err = []
    dt = 0.0001
    T = dt * 4
    N = [4, 8, 12, 16]

    for n in N:
        # Set resolution
        cmd = [
            "mpirun",
            "-np",
            f"{num_processors}",
            "oasismove",
            "NSfracStep",
            f"solver={solver}",
            "problem=TaylorGreen2D",
            "compute_error=1e8",
            f"T={T}",
            f"dt={dt}",
            f"Nx={n}",
            f"Ny={n}",
        ]

        # Run OasisMove
        result = subprocess.run(cmd, capture_output=True, text=True)

        # Assert successful simulation
        assert result.returncode == 0

        # Get output (str)
        output = result.stdout

        match = re.search(
            "Final Error: u0=" + number + " u1=" + number + " p=" + number, str(output)
        )
        err = match.groups()
        u0_err.append(eval(err[0]))
        p_err.append(eval(err[2]))

    p_conv = []
    u_conv = []
    dx = [math.sqrt(2 * (1.0 / n) ** 2) for n in N]
    for i in range(len(dx) - 1):
        u_conv.append(math.log(u0_err[i] / u0_err[i + 1]) / math.log(dx[i] / dx[i + 1]))
        p_conv.append(math.log(p_err[i] / p_err[i + 1]) / math.log(dx[i] / dx[i + 1]))

    assert round(u_conv[-1], 1) == 2.0
    assert round(p_conv[-1], 1) == 4.0


@pytest.mark.parametrize("solver", ["IPCS_ABCN", "IPCS_ABE", "Chorin", "BDFPC_Fast"])
@pytest.mark.parametrize("num_processors", [1, 2])
def test_TaylorGreen2D(num_processors, solver):
    cmd = [
        "mpirun",
        "-np",
        f"{num_processors}",
        "oasismove",
        "NSfracStep",
        f"solver={solver}",
        "problem=TaylorGreen2D",
        "T=0.01",
        "Nx=30",
        "Ny=30",
    ]

    if num_processors > 1 and solver == "Chorin":  # Uses direct solver
        return

    # Run OasisMove
    result = subprocess.run(cmd, capture_output=True, text=True)

    # Assert successful simulation
    assert result.returncode == 0
    output = result.stdout

    match = re.search(
        "Final Error: u0=" + number + " u1=" + number + " p=" + number, str(output)
    )
    err = match.groups()
    if "IPCS_AB" in solver:
        for e in err:
            assert eval(e) < 1e-4
    elif "Chorin" in solver:
        for e in err[:2]:
            assert eval(e) < 1e-4
    else:
        for e in err[:2]:
            assert eval(e) < 1e-4
    # Make sure the optimized version gives the same result as naive
    slow = None
    if "IPCS_AB" in solver:
        slow = "IPCS"
    elif "BDFPC" in solver:
        slow = "BDFPC"
    if slow is not None:
        cmd_slow = [
            "mpirun",
            "-np",
            "1",
            "oasismove",
            "NSfracStep",
            f"solver={slow}",
            "problem=TaylorGreen2D",
            "T=0.01",
            "Nx=30",
            "Ny=30",
        ]

        # Run OasisMove
        result_slow = subprocess.run(cmd_slow, capture_output=True, text=True)

        # Assert successful simulation
        assert result_slow.returncode == 0
        output = result_slow.stdout

        match2 = re.search(
            "Final Error: u0=" + number + " u1=" + number + " p=" + number, str(output)
        )
        err2 = match2.groups()

    if "BDFPC" in solver:
        for e1, e2 in zip(err[:2], err2[:2]):
            assert abs(eval(e1) - eval(e2)) < 1e-9
    elif "IPCS_AB" in solver:
        for e1, e2 in zip(err[:2], err2[:2]):
            assert abs(eval(e1) - eval(e2)) < 1e-8


@pytest.mark.parametrize("solver", ["IPCS_ABCN", "IPCS_ABE"])
@pytest.mark.parametrize("num_processors", [1, 2])
def test_DrivenCavity(num_processors, solver):
    cmd = [
        "mpirun",
        "-np",
        f"{num_processors}",
        "oasismove",
        "NSfracStep",
        "problem=DrivenCavity",
        "T=0.01",
        "Nx=10",
        "Ny=10",
        "plot_interval=10000",
        f"solver={solver}",
        "testing=True",
    ]

    # Run OasisMove
    result = subprocess.run(cmd, capture_output=True, text=True)

    # Assert successful simulation
    assert result.returncode == 0

    # Get output (str)
    output = result.stdout

    match = re.search("Velocity norm = " + number, str(output))
    err = match.groups()

    # Make sure the optimized version gives the same result as naive
    cmd_naive = [
        "mpirun",
        "-np",
        "1",
        "oasismove",
        "NSfracStep",
        "problem=DrivenCavity",
        "T=0.01",
        "Nx=10",
        "Ny=10",
        "plot_interval=10000",
        "solver=IPCS",
        "testing=True",
    ]
    # Run OasisMove
    result_naive = subprocess.run(cmd_naive, capture_output=True, text=True)

    # Assert successful simulation
    assert result_naive.returncode == 0

    # Get output (str)
    output = result_naive.stdout

    match2 = re.search("Velocity norm = " + number, str(output))
    err2 = match2.groups()
    tol = 1e-9 if solver == "IPCS_ABCN" else 5e-3
    assert abs(eval(err[0]) - eval(err2[0])) < tol


if __name__ == "__main__":
    test_DrivenCavity(1)
    test_TaylorGreen2D(1, "IPCS_ABCN")
    test_spatial_rate_of_convergence(1, "BDFPC_Fast")
    test_spatial_rate_of_convergence(1, "IPCS_ABCN")
