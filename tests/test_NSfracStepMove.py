import re
import subprocess

import pytest

number = "([0-9]+.[0-9]+e[+-][0-9]+)"


@pytest.mark.parametrize("num_processors", [1, 2])
def test_DrivenCavity_with_NSfracStepMove(num_processors):
    cmd = ("mpirun -np {} oasism NSfracStepMove problem=DrivenCavity T=0.01 " +
           "Nx=20 Ny=20 plot_interval=10000 solver={} testing=True")
    d = subprocess.check_output(cmd.format(num_processors, "IPCS_ABCN_Move"), shell=True)
    match = re.search("Velocity norm = " + number, str(d))
    err = match.groups()

    velocity_norm = eval(err[0])
    velocity_norm_exact = 4.452349

    assert abs(velocity_norm - velocity_norm_exact) < 1e-16


@pytest.mark.parametrize("num_processors", [1, 2])
def test_MovingVortex(num_processors):
    cmd = ("mpirun -np {} oasism NSfracStepMove problem=MovingVortex T=1 " +
           "Nx=20 Ny=20 solver={}")
    d = subprocess.check_output(cmd.format(num_processors, "IPCS_ABCN_Move"), shell=True)
    match = re.search("Final Error: u0=" + number, str(d))
    err = match.groups()
    assert eval(err[0]) < 5e-3


@pytest.mark.parametrize("num_processors_initial, num_processors_restart", [(1, 2), (2, 2)])
def test_restart_MovingVortex(num_processors_initial, num_processors_restart):
    cmd = ("mpirun -np {} oasism NSfracStepMove problem=MovingVortex T=1 Nx=20 Ny=20 checkpoint=20;" +
           "mpirun -np {} oasism NSfracStepMove problem=MovingVortex T=2 Nx=20 Ny=20 " +
           "restart_folder=results_moving_vortex/data/1/Checkpoint"
           )
    d = subprocess.check_output(cmd.format(num_processors_initial, num_processors_restart), shell=True)
    match = re.findall("u0=" + number, str(d))

    # Check errors from first and second (restart) simulation
    tol_initial = 5E-3
    tol_restart = 5E-4

    assert eval(match[0]) < tol_initial
    assert eval(match[-1]) < tol_restart


if __name__ == '__main__':
    test_MovingVortex(1)
    test_DrivenCavity_with_NSfracStepMove(1)
    test_restart_MovingVortex(1, 1)
