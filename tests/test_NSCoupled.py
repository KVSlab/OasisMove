import re
import subprocess

import pytest

number = "[+-]?([0-9]+.[0-9]+e[+-][0-9]+)"


@pytest.mark.skip(reason="Deprecated test of NSCoupled")
@pytest.mark.parametrize("num_processors", [1])
def test_default_Coupled(num_processors):
    cmd = ["mpirun", "-np", f"{num_processors}", "oasismove", "NSCoupled", "problem=DrivenCavity", "testing=True"]

    # Run OasisMove
    result = subprocess.run(cmd, capture_output=True, text=True)

    # Assert successful simulation
    assert result.returncode == 0

    # Get output (str)
    output = result.stdout

    match = re.search("Velocity in corner = " + number, str(output))
    err = match.groups()
    tol = 1e-16

    assert eval(err[0]) < tol


@pytest.mark.skip(reason="Deprecated test of NSCoupled")
@pytest.mark.parametrize("num_processors", [1])
def test_default_CR_Coupled(num_processors):
    cmd = ["mpirun", "-np", f"{num_processors}", "oasismove", "NSCoupled", "problem=DrivenCavity", "testing=True",
           "element=CR"]

    # Run OasisMove
    result = subprocess.run(cmd, capture_output=True, text=True)

    # Assert successful simulation
    assert result.returncode == 0

    # Get output (str)
    output = result.stdout

    match = re.search("Velocity in corner = " + number, str(output))
    err = match.groups()
    tol = 0.4

    assert eval(err[0]) < tol


if __name__ == '__main__':
    test_default_Coupled()
    test_default_CR_Coupled()
