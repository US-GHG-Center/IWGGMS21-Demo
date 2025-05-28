# syntax=docker/dockerfile:1.4
FROM quay.io/jupyter/julia-notebook:2025-05-12

ENV JULIA_NUM_THREADS=3
ENV OMP_NUM_THREADS=1
ENV XRTM_PATH=/opt/xrtm
# NOTEBOOK_ARGS set to this path will cause the Jupyter Lab environment 
# to start in this folder, so users will see only this Notebook and thus
# know which one to use.
ENV NOTEBOOK_ARGS=/opt/IWGGMS21-Demo/IWGGMS-demonstration.ipynb

USER ${NB_USER}
WORKDIR ${HOME}

# Update the conda enviroment first because the subsequent builds require the
# dependencies.
COPY environment.yml ./
RUN conda env update --name base --quiet --file environment.yml && rm environment.yml && conda clean --all -f -y

USER root
WORKDIR /opt

# Download and install JuliaMono font for good rendering of plotting text
RUN <<EOT
    #!/bin/bash
    set -e
    fonts_dir=/usr/share/fonts/truetype/julia
    curl -OLs https://github.com/cormullion/juliamono/releases/download/v0.060/JuliaMono-ttf.tar.gz
    mkdir -p "${fonts_dir}"
    tar -xzf JuliaMono-ttf.tar.gz -C "${fonts_dir}"
    rm JuliaMono-ttf.tar.gz
EOT

# Clone and make bindx2
RUN <<EOT
    #!/bin/bash
    set -e

    git clone --depth 1 https://github.com/gmcgarragh/bindx2 bindx2
    cp bindx2/make.inc{.example,}
    make -C bindx2
    chown --recursive ${NB_UID}:${NB_GID} bindx2
EOT

# Clone xrtm
RUN <<EOT
    #!/bin/bash
    set -e

    git clone --depth 1 https://github.com/PeterSomkuti/xrtm.git xrtm
    cp xrtm/make.inc{.example,}
EOT

# Patch and make xrtm
COPY make_xrtm.patch ./
RUN <<EOT
    #!/bin/bash
    set -e

    patch -d xrtm <make_xrtm.patch
    rm make_xrtm.patch

    # Various xrtm Makefiles use `${HOME}/bindx2` and `~/bindx2`, so we have to
    # "fudge" HOME to be `/opt`, since that's where we cloned the repo earlier.
    # Making the doc from tex will fail, because, again HOME is used, but we have
    # not cloned tth_C-4.03, so we --ignore-errors to avoid failing.
    HOME=/opt make -C xrtm --ignore-errors

    chown --recursive ${NB_UID}:${NB_GID} xrtm
EOT

# Add the demo notebook into /opt/IWGGMS21-Demo
RUN git clone --depth 1 https://github.com/US-GHG-Center/IWGGMS21-Demo /opt/IWGGMS21-Demo && mkdir /opt/IWGGMS21-Demo/data && chown -R ${NB_USER} /opt/IWGGMS21-Demo


USER ${NB_USER}
WORKDIR ${HOME}

# Install Julia Packages
RUN <<EOT
    julia -e '
        using Pkg;
        Pkg.rm("Pluto"); # Do not need Pluto here (yet)
        Pkg.add(
            [
                "ArgParse",
                "AWS",
                "AWSS3",
                "Dates",
                "DocStringExtensions",
                "HDF5",
                "Interpolations",
                "LaTeXStrings",
                "LinearAlgebra",
                "Logging",
                "LoopVectorization",
                "Plots",
                "Printf",
                "ProgressMeter",
                "Statistics",
                "StatsBase",
                "Unitful"
            ]
        );
        Pkg.add(url="https://github.com/US-GHG-Center/RetrievalToolbox.jl");
        Pkg.precompile();
        Pkg.build();
    '
EOT

# Run the precompile/build within Jupyter (somehow that is necessary)
COPY precompile_build.ipynb ./
RUN jupyter execute precompile_build.ipynb && rm precompile_build.ipynb
