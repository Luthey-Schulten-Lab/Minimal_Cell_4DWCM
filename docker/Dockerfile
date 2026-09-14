###############################################################################
# Minimal Cell 4DWCM — public CUDA image
#
# Build from the repository root:
#   docker build -f docker/Dockerfile -t 4dwcm:ampere .
# Or:  ./docker/build_ampere.sh | ./docker/build_multi.sh | ./docker/build_blackwell.sh
#
# Arch build-args (defaults = Ampere sm_80):
#   CUDA_ARCHITECTURES   LM + LAMMPS CMAKE_CUDA_ARCHITECTURES
#   GPU_ARCH             LAMMPS -DGPU_ARCH (sm_XX)
#   KOKKOS_ARCH_FLAGS    Kokkos -DKokkos_ARCH_*=yes flags
#   HOST_ARCH            Kokkos host arch (default ZEN3)
#   BTREE_NVCC_ARCH      nvcc -arch for btree fused CG
###############################################################################

FROM nvcr.io/nvidia/cuda:12.8.0-devel-ubuntu22.04

ARG CUDA_ARCHITECTURES=80
ARG GPU_ARCH=sm_80
ARG KOKKOS_ARCH_FLAGS="-DKokkos_ARCH_AMPERE80=yes"
ARG HOST_ARCH=ZEN3
ARG BTREE_NVCC_ARCH=sm_80
ARG LM_REPO=https://github.com/Luthey-Schulten-Lab/Lattice_Microbes.git
# Public Lattice Microbes, pinned. The production runs built LM from a private
# 2.6 tree; against that tree's committed state this differs almost entirely by
# formatting, plus one real divergence: MPD_WORDS_PER_SITE is clamped to 16
# there and is MPD_LATTICE_MAX_OCCUPANCY here. Our builds set that to 16, so the
# two are identical. A build with occupancy above 16 would NOT be equivalent.
ARG LM_REF=bc58524c2d76c1b84313672b72b140f7f40b7819
ARG BTREE_REPO=https://github.com/Luthey-Schulten-Lab/btree_chromo_gpu.git
# The commit the 26.8 h production run built from. `protein_science` currently
# points here, but it is a moving branch, so the run is pinned to the SHA.
ARG BTREE_REF=1725d1bb0247df1090baf19e31b3124289ff8964
ARG ODECELL_REPO=https://github.com/Luthey-Schulten-Lab/odecell.git
ARG SC_CHAIN_REPO=https://github.com/Luthey-Schulten-Lab/sc_chain_generation.git
# Commits present in 4dcell-optimize:cuda128, the image the 26.8 h production
# run (job 2413) used. That job rebuilds only LM and btree_chromo, so these two
# are exactly what it ran against.
ARG ODECELL_REF=b19413f2f8c1a88d1fe619e09ad9beb7953ccb96
ARG SC_CHAIN_REF=501a0679f8ab5d0f61f9cda1522b3efbf1800266
ARG LAMMPS_GIT_URL=https://github.com/lammps/lammps.git
# Pinned, not develop. btree_chromo's fix_brownian_kokkos uses the pre-2025
# Kokkos ArrayTypes API (t_x_array, t_v_array, t_f_array_const). Upstream has
# since moved those behind the mixed-precision refactor (KK_ACC_FLOAT), so
# newer LAMMPS fails to compile the DNA additions. This is the release the
# production images were built from.
ARG LAMMPS_GIT_REF=stable_22Jul2025
ARG N_PROC_MAKE=

ENV DEBIAN_FRONTEND=noninteractive
ENV TZ=UTC
SHELL ["/bin/bash", "-c"]

RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        make wget git curl ca-certificates build-essential \
        libstdc++6 libfmt-dev libssl-dev libssh-dev openssh-client \
        libevent-dev zlib1g-dev manpages-dev tzdata gfortran \
        libsundials-dev \
    && rm -rf /var/lib/apt/lists/*

# ---------- Miniconda + LM Python env ----------
RUN curl -sS https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -o /tmp/miniconda.sh && \
    bash /tmp/miniconda.sh -b -p /opt/conda && \
    rm /tmp/miniconda.sh

ENV PATH="/opt/conda/bin:${PATH}"

RUN conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/main && \
    conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/r

COPY docker/lm_precomp.yml /tmp/lm_precomp.yml

RUN source /opt/conda/etc/profile.d/conda.sh && \
    conda env create -n lm_2.5_dev -f /tmp/lm_precomp.yml && \
    conda clean -afy

# ---------- odecell ----------
# Fetched by SHA rather than cloned by branch: `git clone --branch` accepts only
# branch and tag names, so a commit pin needs an explicit shallow fetch.
RUN mkdir -p /Software/odecell && cd /Software/odecell && \
    git init -q && git remote add origin "${ODECELL_REPO}" && \
    git fetch -q --depth 1 origin "${ODECELL_REF}" && \
    git checkout -q FETCH_HEAD

RUN source /opt/conda/etc/profile.d/conda.sh && \
    conda activate lm_2.5_dev && \
    conda install -y -c conda-forge python-libsbml && \
    pip install sbtab==0.9.71 && \
    conda install -y -c conda-forge pycvodes=0.14.0 pydantic=1.8.2 cobra && \
    conda install -y -c anaconda xlrd==1.2.0 mpi4py && \
    conda clean -afy && \
    cd /Software/odecell && \
    pip install .

# ---------- sc_chain_generation ----------
RUN mkdir -p /Software/sc_chain_generation && cd /Software/sc_chain_generation && \
    git init -q && git remote add origin "${SC_CHAIN_REPO}" && \
    git fetch -q --depth 1 origin "${SC_CHAIN_REF}" && \
    git checkout -q FETCH_HEAD

RUN source /opt/conda/etc/profile.d/conda.sh && conda activate lm_2.5_dev && \
    cd /Software/sc_chain_generation/src && make gen_sc_chain && \
    mkdir -p /Software/sc_chain_generation/fortran && \
    cp -f /Software/sc_chain_generation/src/gen_sc_chain \
          /Software/sc_chain_generation/fortran/gen_sc_chain && \
    chmod +x /Software/sc_chain_generation/fortran/gen_sc_chain

# ---------- GCC 11.5 (Kokkos / LAMMPS) ----------
ENV GCC_src=/Software/GCC
ENV GCC_version=11.5.0

RUN mkdir -p ${GCC_src} && \
    wget -O ${GCC_src}/gcc-${GCC_version}.tar.gz \
        https://gcc.gnu.org/pub/gcc/releases/gcc-${GCC_version}/gcc-${GCC_version}.tar.gz && \
    cd ${GCC_src} && tar -xzf gcc-${GCC_version}.tar.gz && \
    cd ${GCC_src}/gcc-${GCC_version} && \
    ./contrib/download_prerequisites && \
    ./configure \
        --build=x86_64-linux-gnu \
        --host=x86_64-linux-gnu \
        --target=x86_64-linux-gnu \
        --enable-checking=release \
        --disable-multilib \
        --enable-languages=c,c++,fortran && \
    make -j"$(nproc)" && make install && \
    cd / && rm -rf ${GCC_src}/gcc-${GCC_version} ${GCC_src}/gcc-${GCC_version}.tar.gz

ENV PATH="/usr/local/bin:${PATH}"
ENV LD_LIBRARY_PATH="/usr/local/lib64:/usr/local/lib:${LD_LIBRARY_PATH}"

ENV cmake_src=/Software/cmake
ENV cmake_version=3.26.4

RUN mkdir -p ${cmake_src} && \
    wget -O ${cmake_src}/cmake-${cmake_version}.tar.gz \
        https://github.com/Kitware/CMake/releases/download/v${cmake_version}/cmake-${cmake_version}.tar.gz && \
    cd ${cmake_src} && tar -xzf cmake-${cmake_version}.tar.gz && \
    cd ${cmake_src}/cmake-${cmake_version} && \
    ./bootstrap && make -j"$(nproc)" && make install && \
    cd / && rm -rf ${cmake_src}/cmake-${cmake_version} ${cmake_src}/cmake-${cmake_version}.tar.gz

ENV OpenMPI_src=/Software/OpenMPI
ENV OpenMPI_version=4.1.5

RUN mkdir -p ${OpenMPI_src} && \
    wget -O ${OpenMPI_src}/openmpi-${OpenMPI_version}.tar.gz \
        https://download.open-mpi.org/release/open-mpi/v4.1/openmpi-${OpenMPI_version}.tar.gz && \
    cd ${OpenMPI_src} && tar -xzf openmpi-${OpenMPI_version}.tar.gz && \
    cd ${OpenMPI_src}/openmpi-${OpenMPI_version} && \
    ./configure \
        --prefix=/usr/local \
        --enable-mpi-cxx \
        CC=gcc CXX=g++ FC=gfortran \
        --with-libevent=internal \
        --enable-static \
        --enable-shared && \
    make -j"$(nproc)" && make install && \
    cd / && rm -rf ${OpenMPI_src}/openmpi-${OpenMPI_version} ${OpenMPI_src}/openmpi-${OpenMPI_version}.tar.gz

# ---------- LAMMPS + DNA package additions from btree_chromo_gpu ----------
ENV LAMMPS_src=/Software/LAMMPS
ENV LAMMPS_tag=OMP_GPU_Kokkos

ENV CC=/usr/local/bin/gcc \
    CXX=/usr/local/bin/g++ \
    MPICC=/usr/local/bin/mpicc \
    MPICXX=/usr/local/bin/mpicxx \
    CUDA_HOME=/usr/local/cuda \
    CUDA_LIB_PATH=/usr/local/cuda/lib64

RUN mkdir -p ${LAMMPS_src} && \
    git clone --depth 1 --branch "${LAMMPS_GIT_REF}" "${LAMMPS_GIT_URL}" "${LAMMPS_src}/lammps"

RUN mkdir -p /Software/btree_chromo && cd /Software/btree_chromo && \
    git init -q && git remote add origin "${BTREE_REPO}" && \
    git fetch -q --depth 1 origin "${BTREE_REF}" && \
    git checkout -q FETCH_HEAD

RUN set -eux; \
    L="${LAMMPS_src}/lammps"; \
    A=/Software/btree_chromo/LAMMPS_src_additions; \
    if [ -d "${A}" ]; then \
      cp -v "${A}/angle_polytorsion.cpp" "${A}/angle_polytorsion.h" \
            "${A}/angle_polytorsionend.cpp" "${A}/angle_polytorsionend.h" \
            "${L}/src/" || true; \
      mkdir -p "${L}/src/OPENMP"; \
      cp -v "${A}/angle_polytorsion_omp.cpp" "${A}/angle_polytorsion_omp.h" \
            "${A}/fix_brownian_base_omp.cpp" "${A}/fix_brownian_base_omp.h" \
            "${A}/fix_brownian_asphere_omp.cpp" "${A}/fix_brownian_asphere_omp.h" \
            "${L}/src/OPENMP/" || true; \
      mkdir -p "${L}/cmake/Modules/Packages"; \
      [ -f "${A}/OPENMP.cmake" ] && cp -v "${A}/OPENMP.cmake" "${L}/cmake/Modules/Packages/OPENMP.cmake" || true; \
      mkdir -p "${L}/src/KOKKOS"; \
      cp -v "${A}/fix_brownian_kokkos.cpp" "${A}/fix_brownian_kokkos.h" \
            "${L}/src/KOKKOS/" || true; \
    fi; \
    if [ -f ${L}/cmake/CMakeLists.txt ]; then \
      sed -i '/# silence nvcc warnings/,/endif()/d' ${L}/cmake/CMakeLists.txt || true; \
    fi; \
    if [ -f ${L}/src/angle_polytorsion.cpp ]; then \
      sed -i \
        's/domain->minimum_image( *\([^,]*\), *\([^,]*\), *\([^)]*\) *);/domain->minimum_image("angle_polytorsion.cpp", __LINE__, \1, \2, \3);/' \
        ${L}/src/angle_polytorsion.cpp || true; \
    fi; \
    if [ -f ${L}/src/angle_polytorsionend.cpp ]; then \
      sed -i \
        's/domain->minimum_image( *\([^,]*\), *\([^,]*\), *\([^)]*\) *);/domain->minimum_image("angle_polytorsionend.cpp", __LINE__, \1, \2, \3);/' \
        ${L}/src/angle_polytorsionend.cpp || true; \
    fi

ENV CMAKE_CUDA_COMPILER=${LAMMPS_src}/lammps/lib/kokkos/bin/nvcc_wrapper

RUN set -eux; \
    NPROC="${N_PROC_MAKE:-$(nproc)}"; \
    cd ${LAMMPS_src} && mkdir -p build_${LAMMPS_tag} && \
    cmake -B build_${LAMMPS_tag} \
          -S ${LAMMPS_src}/lammps/cmake \
          -DCMAKE_BUILD_TYPE=Release \
          -DCMAKE_INSTALL_PREFIX=/usr/local \
          -DBUILD_MPI=yes -DBUILD_OMP=yes -DBUILD_SHARED_LIBS=yes \
          -DLAMMPS_MACHINE=${LAMMPS_tag} \
          -DCMAKE_EXE_LINKER_FLAGS="-L/usr/local/cuda/lib64 -lcuda -lcudart" \
          -DPKG_GPU=yes -DPKG_KOKKOS=on \
          -DKokkos_ENABLE_OMP=yes -DKokkos_ENABLE_CUDA=yes \
          ${KOKKOS_ARCH_FLAGS} \
          -DKokkos_DISABLE_GPU_AUTO_DETECT=on \
          -DKokkos_ARCH_HOSTARCH=${HOST_ARCH} \
          -DGPU_API=cuda -DGPU_ARCH=${GPU_ARCH} \
          -DCMAKE_CUDA_ARCHITECTURES=${CUDA_ARCHITECTURES} \
          -DPKG_OPENMP=yes -DPKG_MOLECULE=yes \
          -DPKG_ASPHERE=yes -DPKG_BROWNIAN=yes \
          -DPKG_RIGID=yes -DPKG_CG-DNA=yes \
          -DPKG_EXTRA-MOLECULE=yes -DPKG_MOLFILE=yes && \
    cd build_${LAMMPS_tag} && make -j "${NPROC}" && make install && \
    rm -rf build_${LAMMPS_tag}

# ---------- btree_chromo (protein_science) ----------
RUN set -eux; \
    NPROC="${N_PROC_MAKE:-$(nproc)}"; \
    cd /Software/btree_chromo && \
    sed -i "s/-arch=sm_[0-9]\\+/-arch=${BTREE_NVCC_ARCH}/" Makefile && \
    make clean || true && \
    make -j "${NPROC}" all

# ---------- Lattice Microbes (public) ----------
ENV LM_dir=/Software/Lattice_Microbes

RUN mkdir -p ${LM_dir} && cd ${LM_dir} && \
    git init -q && git remote add origin "${LM_REPO}" && \
    git fetch -q --depth 1 origin "${LM_REF}" && \
    git checkout -q FETCH_HEAD

RUN set -eux; \
    NPROC="${N_PROC_MAKE:-$(nproc)}"; \
    source /opt/conda/etc/profile.d/conda.sh && conda activate lm_2.5_dev && \
    export CUDA_HOME=/usr/local/cuda && \
    export LD_LIBRARY_PATH=/usr/local/cuda/lib64:${LD_LIBRARY_PATH} && \
    cd ${LM_dir} && rm -rf build && mkdir -p build && cd build && \
    cmake ../src \
        -DMPD_GLOBAL_T_MATRIX=True \
        -DMPD_GLOBAL_R_MATRIX=True \
        -DCUDA_ARCHITECTURES="${CUDA_ARCHITECTURES}" && \
    make -j "${NPROC}" && make install

# ---------- This 4DWCM repository ----------
COPY . /src/4d/
WORKDIR /src/4d

ENV PATH="/opt/conda/envs/lm_2.5_dev/bin:/Software/sc_chain_generation/src:/Software/btree_chromo/build/apps:${PATH}"
ENV PYTHONPATH="/Software/Lattice_Microbes/src/pylm:/Software/Lattice_Microbes/src/jlm:${PYTHONPATH}"
ENV LD_LIBRARY_PATH="/usr/local/cuda/lib64:/usr/local/lib64:/usr/local/lib:${LD_LIBRARY_PATH}"

COPY docker/entrypoint.sh /usr/local/bin/4dwcm-entrypoint.sh
RUN chmod +x /usr/local/bin/4dwcm-entrypoint.sh

RUN echo "source /opt/conda/etc/profile.d/conda.sh" >> /root/.bashrc && \
    echo "conda activate lm_2.5_dev" >> /root/.bashrc

ENTRYPOINT ["/usr/local/bin/4dwcm-entrypoint.sh"]
CMD ["bash"]
