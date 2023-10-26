import re
import subprocess
from os import getcwd, path

import pytest

number = "([0-9]+.[0-9]+e[+-][0-9]+)"
number_pattern = r"(\d+.\d+)"


@pytest.mark.parametrize("num_processors", [1, 2])
def test_DrivenCavity_with_NSfracStepMove(num_processors):
    cmd = ["mpirun", "-np", f"{num_processors}", "oasismove", "NSfracStepMove", "problem=DrivenCavity", "T=0.01",
           "Nx=10", "Ny=10", "plot_interval=10000", "solver=IPCS_ABCN_Move", "testing=True"]

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


@pytest.mark.parametrize("num_processors", [1])
def test_MovingCylinder(num_processors):
    # Simulation parameters
    dt = 0.001
    T = 10 * dt
    cwd = getcwd()
    mesh_path = path.join(cwd, "src/oasismove/mesh/cylinder.xdmf")

    cmd = ["mpirun", "-np", f"{num_processors}", "oasismove", "NSfracStepMove", "problem=MovingCylinder",
           f"T={T}", f"dt={dt}", f"mesh_path={mesh_path}", "checkpoint=10", "solver=IPCS_ABCN_Move", "save_step=10"]

    # Run OasisMove
    result = subprocess.run(cmd, capture_output=True, text=True)

    # Assert successful simulation
    assert result.returncode == 0

    # Get output (str)
    output = result.stdout

    pattern = re.compile(r"velocity=" + number_pattern)
    velocities = []
    for match in pattern.finditer(str(output)):
        velocities.append(eval(match.group(1)))

    u_exact_max = 2.016
    u_exact_mean = 1.004
    tol = 1E-12

    assert abs(u_exact_max - velocities[0]) < tol
    assert abs(u_exact_mean - velocities[1]) < tol


@pytest.mark.parametrize("num_processors", [1])
def test_MovingWall(num_processors):
    # Simulation parameters
    dt = 0.01
    T = 10 * dt

    cmd = ["mpirun", "-np", f"{num_processors}", "oasismove", "NSfracStepMove", "problem=MovingWall", f"T={T}",
           f"dt={dt}", "Nx=10", "Ny=10", "solver=IPCS_ABCN_Move"]

    # Run OasisMove
    result = subprocess.run(cmd, capture_output=True, text=True)

    # Assert successful simulation
    assert result.returncode == 0
    # Get output (str)
    output = result.stdout

    pattern = re.compile(r"Reynolds number=" + number_pattern)
    reynolds_numbers = []
    for match in pattern.finditer(str(output)):
        reynolds_numbers.append(eval(match.group(1)))

    re_exact_max = 6.698
    re_exact_mean = 2.787
    tol = 1E-12

    assert abs(re_exact_max - reynolds_numbers[0]) < tol
    assert abs(re_exact_mean - reynolds_numbers[1]) < tol


@pytest.mark.parametrize("num_processors", [1])
def test_MovingTaylorGreen3D(num_processors):
    # Simulation parameters
    dt = 0.001
    T = 10 * dt

    cmd = ["mpirun", "-np", f"{num_processors}", "oasismove", "NSfracStepMove", "problem=MovingTaylorGreen3D", f"T={T}",
           f"dt={dt}", "Nx=10", "Ny=10", "Nz=10", "solver=IPCS_ABCN_Move"]

    # Run OasisMove
    result = subprocess.run(cmd, capture_output=True, text=True)

    # Assert successful simulation
    assert result.returncode == 0
    # Get output (str)
    output = result.stdout

    pattern = re.compile(r"Reynolds number=" + number_pattern)
    reynolds_numbers = []
    for match in pattern.finditer(str(output)):
        reynolds_numbers.append(eval(match.group(1)))

    re_exact_max = 8669.186
    re_exact_mean = 3964.611
    tol = 1E-12

    assert abs(re_exact_max - reynolds_numbers[0]) < tol
    assert abs(re_exact_mean - reynolds_numbers[1]) < tol


if __name__ == '__main__':
    test_MovingVortex(1)
    test_MovingVortex(2)
    test_DrivenCavity_with_NSfracStepMove(1)
    test_DrivenCavity_with_NSfracStepMove(2)
    test_MovingWall(1)
    test_MovingCylinder(1)
    test_MovingTaylorGreen3D(1)
