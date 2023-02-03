# Installation
## Installing OasisMove using `conda`

We recommend installing *OasisMove* and its dependencies through `conda`.  
Start by cloning into the repository:

``` console
$ git clone https://github.com/KVSLab/OasisMove.git
$ cd OasisMove
```

Then, using the ``environment.yml`` file in the root of the repository, you can call:

``` console
$ conda env update --file environment.yml --name your_environment
```

Next, can now activate your environment by running::

``` console
$ conda activate your_environment
```

Finally, install OasisMove inside your `conda` environment using `pip`: 

``` console
$ python3 -m pip install .
```

Now you are all set, and can start using OasisMove.

## Editable installation

If you want to make changes to the underlying OasisMove package, you can install an editable version by adding the `--editable` flag:

``` console
$ python3 -m pip install --editable .
```

The ``--editable`` flag installs the project in editable mode meaning that any changes to the original package will be
reflected directly in your environment.

## Running OasisMove using Docker

A `Dockerfile` is supplied in the root directory of the repository, which can build a docker-image with all
dependencies installed. The Docker-image can be built with the following command:

``` console
$ docker build -t name_of_image -f docker/Dockerfile .
```

A Docker-container can then be started with the following command:

``` console
$ docker run -ti --network=host -e DISPLAY=$DISPLAY -v /tmp/.X11-unix:/tmp/.X11-unix -v $(pwd):/root/shared -w /root/shared --rm --shm-size=512m name_of_image
```

To run the OasisMove GUI, you need to call:

``` console
$ xhost +local:root
```

on your system before running the scripts.

### Note on Docker

Remember to call:

``` console
xhost -local:root
```

on the host system when you are done running the Docker container.
    