#!/usr/bin/env python

import os
import sys

sys.path.append(os.getcwd())


def main():
    assert sys.argv[1] in ('NSfracStep', 'NSCoupled', 'NSfracStepMove')
    solver = sys.argv.pop(1)
    if solver == 'NSfracStep':
        from oasismove import NSfracStep

    elif solver == 'NSCoupled':
        from oasismove import NSCoupled

    elif solver == "NSfracStepMove":
        from oasismove import NSfracStepMove

    else:
        print(sys.argv[1])
        raise NotImplementedError


if __name__ == '__main__':
    main()
