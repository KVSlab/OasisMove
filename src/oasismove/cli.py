"""Console script for oasismove."""
import argparse


def main():
    """Console script for oasismove."""
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("_", nargs="*")
    args = parser.parse_args()

    print("Arguments: " + str(args._))
    print("Replace this message by putting your code into oasismove.cli.main")
    return 0
