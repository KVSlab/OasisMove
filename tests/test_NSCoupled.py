import re
import subprocess

number = "[+-]?([0-9]+.[0-9]+e[+-][0-9]+)"


def test_default_Coupled():
    # Set resolution
    Nx = Ny = 50
    is_testing = True

    # Command
    cmd = (f"mpirun -np 1 oasism NSCoupled problem=DrivenCavity Nx={Nx} Ny={Ny} testing={is_testing}")

    d = subprocess.check_output(cmd, shell=True)

    match = re.search("Velocity in corner = " + number, str(d))
    err = match.groups()
    err_value = eval(err[0])
    tol = 1e-16

    assert abs(err_value) < tol


if __name__ == '__main__':
    test_default_Coupled()
