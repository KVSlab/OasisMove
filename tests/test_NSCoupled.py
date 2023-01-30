import re
import subprocess

number = "[+-]?([0-9]+.[0-9]+e[+-][0-9]+)"


def test_default_Coupled():
    d = subprocess.check_output("mpirun -np 1 oasism NSCoupled problem=DrivenCavity testing=True", shell=True)
    match = re.search("Velocity in corner = " + number, str(d))
    err = match.groups()
    tol = 1e-16

    assert eval(err[0]) < tol


def test_default_CR_Coupled():
    d = subprocess.check_output("mpirun -np 1 oasism NSCoupled problem=DrivenCavity testing=True element=CR",
                                shell=True)
    match = re.search("Velocity in corner = " + number, str(d))
    err = match.groups()
    tol = 0.4

    assert eval(err[0]) < tol


if __name__ == '__main__':
    test_default_Coupled()
    test_default_CR_Coupled()
