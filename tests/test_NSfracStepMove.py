import re
import subprocess

import pytest

number = "([0-9]+.[0-9]+e[+-][0-9]+)"


@pytest.mark.parametrize("num_processors", [1, 2])
def test_DrivenCavity_with_NSfracStepMove(num_processors):
    cmd = ["mpirun", "-np", f"{num_processors}", "oasismove", "NSfracStepMove", "problem=DrivenCavity", "T=0.01"
        , "Nx=10", "Ny=10", "plot_interval=10000", "solver=IPCS_ABCN_Move", "testing=True"]

    # Run OasisMove
    result = subprocess.run(cmd, capture_output=True, text=True)

    # Assert successful simulation
    assert result.returncode == 0

    # Get output (str)
    output = result.stdout

    match = re.search("Velocity norm = " + number, str(output))
    err = match.groups()

    velocity_norm = eval(err[0])
    velocity_norm_exact = 3.117529

    assert abs(velocity_norm - velocity_norm_exact) < 1e-16


@pytest.mark.parametrize("num_processors", [1, 2])
def test_MovingVortex(num_processors):
    cmd = ["mpirun", "-np", f"{num_processors}", "oasismove", "NSfracStepMove", "problem=MovingVortex", "T=1",
           "Nx=10", "Ny=10", "solver=IPCS_ABCN_Move"]

    # Run OasisMove
    result = subprocess.run(cmd, capture_output=True, text=True)

    # Assert successful simulation
    assert result.returncode == 0

    # Get output (str)
    output = result.stdout

    match = re.search("Final Error: u0=" + number, str(output))
    err = match.groups()
    assert eval(err[0]) < 5e-3


if __name__ == '__main__':
    test_MovingVortex(1)
    test_MovingVortex(2)
    test_DrivenCavity_with_NSfracStepMove(1)
    test_DrivenCavity_with_NSfracStepMove(2)
